#!/usr/bin/env python3
"""Conservatively rescue dual-index reads from demultiplexed Undetermined FASTQ.

The script never edits raw FASTQ files.  Assigned FASTQs are symlinked into a new
directory and rescued reads are written as extra FASTQ chunks named for the
corresponding Sample_ID.  Only a unique nearest expected i7+i5 pair within the
configured mismatch limit is rescued; ties and malformed/missing indexes remain
Undetermined.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from contextlib import ExitStack
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import combinations, product, zip_longest
from pathlib import Path
from typing import Callable, Iterable, Iterator, TextIO


FASTQ_SUFFIX_RE = re.compile(r"\.(?:fastq|fq)(?:\.gz)?$", re.IGNORECASE)
READ_TOKEN_RE = re.compile(
    r"^(?P<prefix>.+?)(?:_S\d+)?(?:_L(?P<lane>\d{3}))?"
    r"[_\.](?P<read>R?[12])(?:[_\.](?P<chunk>\d+))?$",
    re.IGNORECASE,
)
DNA_RE = re.compile(r"^[ACGTN]+$", re.IGNORECASE)
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
DNA_SYMBOLS = "ACGTN"


@dataclass(frozen=True)
class FastqRecord:
    header: str
    sequence: str
    plus: str
    quality: str


@dataclass(frozen=True)
class FastqFile:
    path: Path
    sample: str
    read: str
    lane: str
    chunk: str
    key: str


@dataclass(frozen=True)
class ExpectedIndex:
    sample: str
    i7: str
    i5: str


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1].upper()


def open_text(path: Path, mode: str = "rt", compresslevel: int = 1):
    if path.suffix.lower() == ".gz":
        return gzip.open(
            path,
            mode,
            encoding=None if "b" in mode else "ascii",
            compresslevel=compresslevel,
        )
    return path.open(mode, encoding=None if "b" in mode else "ascii")


def iter_fastq(
    path: Path,
    compressed_position_callback: Callable[[int], None] | None = None,
) -> Iterator[FastqRecord]:
    raw_handle = None
    stack = ExitStack()
    try:
        if path.suffix.lower() == ".gz" and compressed_position_callback is not None:
            raw_handle = stack.enter_context(path.open("rb"))
            gzip_handle = stack.enter_context(gzip.GzipFile(fileobj=raw_handle, mode="rb"))
            handle = stack.enter_context(io.TextIOWrapper(gzip_handle, encoding="ascii"))
        else:
            handle = stack.enter_context(open_text(path, "rt"))
        line_number = 0
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            line_number += 4
            if not sequence or not plus or not quality:
                raise ValueError(f"truncated FASTQ near line {line_number}: {path}")
            record = FastqRecord(
                header.rstrip("\r\n"),
                sequence.rstrip("\r\n").upper(),
                plus.rstrip("\r\n"),
                quality.rstrip("\r\n"),
            )
            if not record.header.startswith("@") or not record.plus.startswith("+"):
                raise ValueError(f"invalid FASTQ record near line {line_number}: {path}")
            if len(record.sequence) != len(record.quality):
                raise ValueError(f"sequence/quality length mismatch near line {line_number}: {path}")
            if compressed_position_callback is not None:
                position = raw_handle.tell() if raw_handle is not None else handle.tell()
                compressed_position_callback(position)
            yield record
    finally:
        stack.close()


def write_record(handle: TextIO, record: FastqRecord) -> None:
    handle.write(record.header + "\n")
    handle.write(record.sequence + "\n")
    handle.write(record.plus + "\n")
    handle.write(record.quality + "\n")


def canonical_read_id(header: str) -> str:
    return re.sub(r"/[12]$", "", header.split()[0])


def extract_header_indexes(header: str) -> tuple[str, str] | None:
    """Extract i7+i5 from a CASAVA/BCL Convert FASTQ header."""
    parts = header.split()
    if len(parts) < 2:
        return None
    candidate = parts[-1].split(":")[-1].upper().replace("-", "+").replace("_", "+")
    pieces = candidate.split("+")
    if len(pieces) != 2:
        return None
    i7, i5 = pieces
    if not (DNA_RE.fullmatch(i7) and DNA_RE.fullmatch(i5)):
        return None
    return i7, i5


def parse_fastq_name(path: Path) -> FastqFile | None:
    stem = FASTQ_SUFFIX_RE.sub("", path.name)
    match = READ_TOKEN_RE.match(stem)
    if not match:
        return None
    prefix = match.group("prefix")
    read = match.group("read").upper()
    if read in {"1", "2"}:
        read = "R" + read
    lane = match.group("lane") or ""
    chunk = match.group("chunk") or ""
    sample = re.sub(r"_S\d+$", "", prefix)
    return FastqFile(path, sample, read, lane, chunk, f"{sample}|{lane}|{chunk}")


def discover_fastqs(input_dir: Path) -> tuple[list[FastqFile], list[str]]:
    parsed: list[FastqFile] = []
    warnings: list[str] = []
    for path in sorted(input_dir.rglob("*")):
        if not path.is_file() or not FASTQ_SUFFIX_RE.search(path.name):
            continue
        item = parse_fastq_name(path)
        if item is None:
            warnings.append(f"Unrecognized FASTQ filename: {path}")
        else:
            parsed.append(item)
    return parsed, warnings


def read_sample_sheet(path: Path) -> list[ExpectedIndex]:
    lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    data_start = 0
    sectioned = False
    data_sections = {"[data]", "[bclconvert_data]"}
    for number, line in enumerate(lines):
        if line.strip().lower() in data_sections:
            data_start = number + 1
            sectioned = True
            break
    data_lines: list[str] = []
    for line in lines[data_start:]:
        stripped = line.strip()
        if not stripped:
            continue
        if sectioned and stripped.startswith("[") and stripped.endswith("]"):
            break
        data_lines.append(line)
    if not data_lines:
        raise ValueError(f"SampleSheet contains no [Data] rows: {path}")
    reader = csv.DictReader(data_lines)
    if reader.fieldnames is None:
        raise ValueError(f"SampleSheet has no header: {path}")
    normalized = {name.strip().lower(): name for name in reader.fieldnames}
    sample_column = next(
        (normalized[key] for key in ("sample_id", "sampleid", "sample_name") if key in normalized),
        None,
    )
    i7_column = next(
        (normalized[key] for key in ("index", "i7", "index1") if key in normalized),
        None,
    )
    i5_column = next(
        (normalized[key] for key in ("index2", "i5") if key in normalized),
        None,
    )
    if not sample_column or not i7_column or not i5_column:
        raise ValueError("SampleSheet must contain Sample_ID, index, and index2")
    records: list[ExpectedIndex] = []
    for row in reader:
        sample = (row.get(sample_column) or "").strip()
        i7 = re.sub(r"[^A-Za-z]", "", row.get(i7_column) or "").upper()
        i5 = re.sub(r"[^A-Za-z]", "", row.get(i5_column) or "").upper()
        if not sample and not i7 and not i5:
            continue
        if not sample or not DNA_RE.fullmatch(i7) or not DNA_RE.fullmatch(i5):
            raise ValueError(f"invalid SampleSheet row for sample={sample!r}")
        records.append(ExpectedIndex(sample, i7, i5))
    if len(records) != 7:
        raise ValueError(f"exactly 7 indexed samples are required; found {len(records)}")
    samples = [record.sample for record in records]
    pairs = [(record.i7, record.i5) for record in records]
    if len(set(samples)) != 7 or len(set(pairs)) != 7:
        raise ValueError("Sample_ID and dual-index pairs must each be unique")
    return records


def hamming(left: str, right: str) -> int | None:
    if len(left) != len(right):
        return None
    return sum(a != b for a, b in zip(left, right))


def orient_expected(records: Iterable[ExpectedIndex], orientation: str) -> list[ExpectedIndex]:
    if orientation not in {"as_given", "reverse_complement"}:
        raise ValueError(f"invalid i5 orientation: {orientation}")
    return [
        ExpectedIndex(
            record.sample,
            record.i7,
            record.i5 if orientation == "as_given" else reverse_complement(record.i5),
        )
        for record in records
    ]


def classify_index(
    observed: tuple[str, str] | None,
    expected: list[ExpectedIndex],
    max_total_mismatches: int,
    max_per_index_mismatches: int,
    minimum_margin: int,
) -> tuple[str | None, str, int | None, int | None, int | None]:
    if observed is None:
        return None, "missing_or_malformed_index", None, None, None
    i7, i5 = observed
    candidates: list[tuple[int, int, int, str]] = []
    for record in expected:
        d7 = hamming(i7, record.i7)
        d5 = hamming(i5, record.i5)
        if d7 is None or d5 is None:
            continue
        candidates.append((d7 + d5, d7, d5, record.sample))
    if not candidates:
        return None, "index_length_mismatch", None, None, None
    candidates.sort()
    best_total, best_d7, best_d5, best_sample = candidates[0]
    if best_total > max_total_mismatches or max(best_d7, best_d5) > max_per_index_mismatches:
        return None, "too_distant", best_total, best_d7, best_d5
    if len(candidates) > 1:
        second_total = candidates[1][0]
        if second_total == best_total or second_total - best_total < minimum_margin:
            return None, "ambiguous_nearest_index", best_total, best_d7, best_d5
    return best_sample, f"rescued_distance_{best_total}", best_total, best_d7, best_d5


def sequence_neighbors(sequence: str, max_mismatches: int) -> Iterator[tuple[str, int]]:
    """Yield all DNA strings within a small Hamming radius of ``sequence``."""
    sequence = sequence.upper()
    yield sequence, 0
    for distance in range(1, max_mismatches + 1):
        for positions in combinations(range(len(sequence)), distance):
            replacements = [
                tuple(symbol for symbol in DNA_SYMBOLS if symbol != sequence[position])
                for position in positions
            ]
            for symbols in product(*replacements):
                value = list(sequence)
                for position, symbol in zip(positions, symbols):
                    value[position] = symbol
                yield "".join(value), distance


class FastIndexClassifier:
    """O(1) classification for the configured rescue radius, with fallback caching.

    A seven-sample experiment with an allowed dual-index distance of two has only
    a few tens of thousands of relevant observed pairs. Precomputing those pairs
    avoids repeating Python Hamming-distance loops for every FASTQ record.
    """

    def __init__(
        self,
        expected: list[ExpectedIndex],
        max_total_mismatches: int,
        max_per_index_mismatches: int,
        minimum_margin: int,
        cache_max_entries: int = 500_000,
    ):
        self.expected = expected
        self.max_total_mismatches = max_total_mismatches
        self.max_per_index_mismatches = max_per_index_mismatches
        self.minimum_margin = minimum_margin
        self.cache_max_entries = cache_max_entries
        self.lookup: dict[
            tuple[str, str], tuple[str | None, str, int | None, int | None, int | None]
        ] = {}
        self.fallback_cache: dict[
            tuple[str, str], tuple[str | None, str, int | None, int | None, int | None]
        ] = {}
        self.counters: Counter[str] = Counter()
        self._build_lookup()

    def _build_lookup(self) -> None:
        # Larger radii grow combinatorially. The normal, validated pipeline uses
        # radius two; unusual settings keep exact semantics via the cached fallback.
        if self.max_total_mismatches > 2:
            return
        candidates: set[tuple[str, str]] = set()
        per_index_radius = min(self.max_total_mismatches, self.max_per_index_mismatches)
        for record in self.expected:
            i7_neighbors = list(sequence_neighbors(record.i7, per_index_radius))
            i5_neighbors = list(sequence_neighbors(record.i5, per_index_radius))
            for (i7, d7), (i5, d5) in product(i7_neighbors, i5_neighbors):
                if d7 + d5 <= self.max_total_mismatches:
                    candidates.add((i7, i5))
        self.lookup = {
            observed: classify_index(
                observed,
                self.expected,
                self.max_total_mismatches,
                self.max_per_index_mismatches,
                self.minimum_margin,
            )
            for observed in candidates
        }

    def classify(
        self, observed: tuple[str, str] | None
    ) -> tuple[str | None, str, int | None, int | None, int | None]:
        self.counters["calls"] += 1
        if observed is None:
            self.counters["missing_index"] += 1
            return classify_index(
                observed,
                self.expected,
                self.max_total_mismatches,
                self.max_per_index_mismatches,
                self.minimum_margin,
            )
        result = self.lookup.get(observed)
        if result is not None:
            self.counters["precomputed_lookup_hits"] += 1
            return result
        result = self.fallback_cache.get(observed)
        if result is not None:
            self.counters["fallback_cache_hits"] += 1
            return result
        self.counters["full_distance_calculations"] += 1
        result = classify_index(
            observed,
            self.expected,
            self.max_total_mismatches,
            self.max_per_index_mismatches,
            self.minimum_margin,
        )
        if len(self.fallback_cache) < self.cache_max_entries:
            self.fallback_cache[observed] = result
        return result

    def stats(self) -> dict[str, int]:
        return {
            "precomputed_lookup_entries": len(self.lookup),
            "fallback_cache_entries": len(self.fallback_cache),
            **dict(self.counters),
        }


def choose_i5_orientation(
    parsed: list[FastqFile],
    records: list[ExpectedIndex],
    scan_reads_per_file: int,
    max_total_mismatches: int,
) -> tuple[str, dict[str, int]]:
    scores = {"as_given": 0, "reverse_complement": 0}
    observed_headers = 0
    r1_files = [item for item in parsed if item.read == "R1"]
    for item in r1_files:
        for number, record in enumerate(iter_fastq(item.path), start=1):
            observed = extract_header_indexes(record.header)
            if observed is not None:
                observed_headers += 1
                for orientation in scores:
                    sample, _, _, _, _ = classify_index(
                        observed,
                        orient_expected(records, orientation),
                        max_total_mismatches,
                        max_total_mismatches,
                        1,
                    )
                    scores[orientation] += int(sample is not None)
            if number >= scan_reads_per_file:
                break
    if observed_headers == 0:
        raise ValueError(
            "No i7+i5 strings were found in FASTQ headers. Read-level rescue is impossible "
            "without synchronized I1/I2 FASTQ or the original BCL run folder."
        )
    selected = "reverse_complement" if scores["reverse_complement"] > scores["as_given"] else "as_given"
    return selected, {**scores, "headers_scanned_with_dual_index": observed_headers}


def pair_undetermined_files(parsed: list[FastqFile]) -> list[tuple[FastqFile, FastqFile | None]]:
    by_key: dict[str, dict[str, FastqFile]] = defaultdict(dict)
    for item in parsed:
        if item.sample.lower().startswith("undetermined"):
            by_key[item.key][item.read] = item
    pairs: list[tuple[FastqFile, FastqFile | None]] = []
    for key in sorted(by_key):
        r1 = by_key[key].get("R1")
        r2 = by_key[key].get("R2")
        if r1 is None:
            raise ValueError(f"Undetermined chunk has R2 without R1: {key}")
        pairs.append((r1, r2))
    return pairs


def link_assigned_fastqs(parsed: list[FastqFile], output_fastq_dir: Path) -> int:
    linked = 0
    output_fastq_dir.mkdir(parents=True, exist_ok=True)
    for item in parsed:
        if item.sample.lower().startswith("undetermined"):
            continue
        destination = output_fastq_dir / item.path.name
        if destination.exists() or destination.is_symlink():
            raise FileExistsError(f"output FASTQ already exists: {destination}")
        destination.symlink_to(item.path.resolve())
        linked += 1
    return linked


class ProcessGzipWriter:
    """Text writer backed by a dedicated gzip/pigz subprocess."""

    def __init__(
        self,
        path: Path,
        backend: str,
        compresslevel: int,
        pigz_threads: int,
    ):
        self.path = path
        self.backend = backend
        self.output_handle = None
        if backend == "pigz":
            pigz = shutil.which("pigz")
            if pigz is None:
                raise RuntimeError("compression backend 'pigz' requested but pigz is not installed")
            self.output_handle = path.open("wb")
            command = [
                pigz,
                "-c",
                "-p",
                str(pigz_threads),
                f"-{compresslevel}",
            ]
            self.process = subprocess.Popen(
                command,
                stdin=subprocess.PIPE,
                stdout=self.output_handle,
            )
        else:
            worker = Path(__file__).with_name("gzip_stream_worker.py")
            command = [
                sys.executable,
                str(worker),
                "--output",
                str(path),
                "--compresslevel",
                str(compresslevel),
            ]
            self.process = subprocess.Popen(command, stdin=subprocess.PIPE)
        if self.process.stdin is None:
            raise RuntimeError(f"failed to open compressor stdin for {path}")
        self.text = io.TextIOWrapper(
            self.process.stdin,
            encoding="ascii",
            newline="",
            write_through=False,
        )

    def close_input(self) -> None:
        if not self.text.closed:
            self.text.close()

    def wait(self) -> None:
        return_code = self.process.wait()
        if self.output_handle is not None:
            self.output_handle.close()
        if return_code != 0:
            raise RuntimeError(
                f"{self.backend} compressor failed with exit code {return_code}: {self.path}"
            )


class LazyFastqWriters:
    def __init__(
        self,
        output_dir: Path,
        chunk_number: int,
        paired: bool,
        compresslevel: int,
        compression_backend: str,
        pigz_threads: int,
    ):
        self.output_dir = output_dir
        self.chunk_number = chunk_number
        self.paired = paired
        self.compresslevel = compresslevel
        if compression_backend == "auto":
            compression_backend = "pigz" if shutil.which("pigz") else "parallel_python"
        if compression_backend == "pigz" and shutil.which("pigz") is None:
            raise RuntimeError("RESCUE_COMPRESSION_BACKEND=pigz but pigz was not found")
        self.compression_backend = compression_backend
        self.pigz_threads = pigz_threads
        self.stack = ExitStack()
        self.handles: dict[tuple[str, str], TextIO] = {}
        self.process_writers: list[ProcessGzipWriter] = []

    def _handle(self, sample: str, read: str) -> TextIO:
        key = (sample, read)
        if key not in self.handles:
            safe_sample = re.sub(r"[^A-Za-z0-9_.-]+", "_", sample).strip("_")
            filename = f"{safe_sample}_L900_{read}_{self.chunk_number:03d}.fastq.gz"
            path = self.output_dir / filename
            if self.compression_backend == "gzip":
                handle = gzip.open(path, "wt", encoding="ascii", compresslevel=self.compresslevel)
                self.handles[key] = self.stack.enter_context(handle)
            else:
                process_writer = ProcessGzipWriter(
                    path,
                    self.compression_backend,
                    self.compresslevel,
                    self.pigz_threads,
                )
                self.process_writers.append(process_writer)
                self.handles[key] = process_writer.text
        return self.handles[key]

    def write(self, sample: str, record1: FastqRecord, record2: FastqRecord | None) -> None:
        write_record(self._handle(sample, "R1"), record1)
        if record2 is not None:
            write_record(self._handle(sample, "R2"), record2)

    def close(self) -> None:
        # Signal EOF to every compressor before waiting. This avoids serial waits
        # and lets all sample/read outputs finish concurrently.
        for writer in self.process_writers:
            writer.close_input()
        for writer in self.process_writers:
            writer.wait()
        self.stack.close()


def format_duration(seconds: float | None) -> str:
    if seconds is None or seconds < 0:
        return "unknown"
    rounded = int(round(seconds))
    hours, remainder = divmod(rounded, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours:
        return f"{hours:d}h {minutes:02d}m {secs:02d}s"
    if minutes:
        return f"{minutes:d}m {secs:02d}s"
    return f"{secs:d}s"


class RescueProgressReporter:
    """Periodically report rescue progress without external dependencies."""

    def __init__(
        self,
        output_path: Path,
        total_compressed_r1_bytes: int,
        interval_seconds: float,
    ):
        self.output_path = output_path
        self.total_bytes = max(total_compressed_r1_bytes, 1)
        self.interval_seconds = interval_seconds
        self.started_monotonic = time.monotonic()
        self.started_at = datetime.now(timezone.utc).isoformat()
        self.last_report_monotonic = 0.0

    def report(
        self,
        processed: int,
        rescued: int,
        bytes_done: int,
        chunk_number: int,
        chunk_count: int,
        force: bool = False,
        status: str = "running",
    ) -> None:
        now = time.monotonic()
        if not force and now - self.last_report_monotonic < self.interval_seconds:
            return
        elapsed = max(now - self.started_monotonic, 1e-9)
        bounded_bytes = min(max(bytes_done, 0), self.total_bytes)
        fraction = bounded_bytes / self.total_bytes
        read_rate = processed / elapsed
        byte_rate = bounded_bytes / elapsed
        eta_seconds = (
            (self.total_bytes - bounded_bytes) / byte_rate
            if byte_rate > 0 and bounded_bytes < self.total_bytes
            else 0.0 if bounded_bytes >= self.total_bytes else None
        )
        rescue_percent = 100 * rescued / processed if processed else 0.0
        payload = {
            "status": status,
            "updated_at": datetime.now(timezone.utc).isoformat(),
            "started_at": self.started_at,
            "chunk": chunk_number,
            "chunks_total": chunk_count,
            "processed_read_pairs_or_reads": processed,
            "rescued_read_pairs_or_reads": rescued,
            "rescued_percent_so_far": rescue_percent,
            "compressed_r1_bytes_processed": bounded_bytes,
            "compressed_r1_bytes_total": self.total_bytes,
            "progress_percent": 100 * fraction,
            "read_pairs_or_reads_per_second": read_rate,
            "elapsed_seconds": elapsed,
            "eta_seconds": eta_seconds,
        }
        temporary = self.output_path.with_suffix(self.output_path.suffix + ".tmp")
        temporary.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        temporary.replace(self.output_path)
        eta_text = format_duration(eta_seconds)
        print(
            f"[rescue] {100 * fraction:6.2f}% | chunk {chunk_number}/{chunk_count} | "
            f"processed {processed:,} | rescued {rescued:,} ({rescue_percent:.2f}%) | "
            f"{read_rate:,.0f} read pairs/s | elapsed {format_duration(elapsed)} | "
            f"ETA {eta_text}",
            file=sys.stderr,
            flush=True,
        )
        self.last_report_monotonic = now


def minimum_expected_pair_distance(records: list[ExpectedIndex]) -> int:
    distances: list[int] = []
    for left_index, left in enumerate(records):
        for right in records[left_index + 1 :]:
            d7 = hamming(left.i7, right.i7)
            d5 = hamming(left.i5, right.i5)
            if d7 is not None and d5 is not None:
                distances.append(d7 + d5)
    return min(distances) if distances else 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Rescue dual-index Undetermined FASTQ reads")
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--sample-sheet", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--mode", choices=["inspect", "rescue"], default="rescue")
    parser.add_argument("--max-total-mismatches", type=int, default=2)
    parser.add_argument("--max-per-index-mismatches", type=int, default=2)
    parser.add_argument("--minimum-margin", type=int, default=1)
    parser.add_argument(
        "--i5-orientation",
        choices=["auto", "as_given", "reverse_complement"],
        default="auto",
    )
    parser.add_argument("--orientation-scan-reads", type=int, default=100_000)
    parser.add_argument("--compresslevel", type=int, choices=range(1, 10), default=1)
    parser.add_argument(
        "--compression-backend",
        choices=["auto", "pigz", "parallel_python", "gzip"],
        default="auto",
        help="Output compression: auto uses pigz when available, otherwise subprocess gzip",
    )
    parser.add_argument(
        "--pigz-threads-per-file",
        type=int,
        default=4,
        help="pigz threads per active output FASTQ (default: 4)",
    )
    parser.add_argument(
        "--progress-interval-seconds",
        type=float,
        default=30.0,
        help="Seconds between rescue progress updates (default: 30)",
    )
    parser.add_argument(
        "--progress-check-reads",
        type=int,
        default=100_000,
        help="Check whether a progress update is due every N reads/read pairs",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.max_total_mismatches < 0 or args.max_per_index_mismatches < 0:
        raise ValueError("mismatch limits cannot be negative")
    if args.minimum_margin < 1:
        raise ValueError("--minimum-margin must be at least 1")
    if args.progress_interval_seconds < 0:
        raise ValueError("--progress-interval-seconds cannot be negative")
    if args.progress_check_reads < 1:
        raise ValueError("--progress-check-reads must be at least 1")
    if args.pigz_threads_per_file < 1:
        raise ValueError("--pigz-threads-per-file must be at least 1")
    parsed, warnings = discover_fastqs(args.input_dir)
    if not parsed:
        raise ValueError(f"no R1/R2 FASTQ files found under {args.input_dir}")
    records = read_sample_sheet(args.sample_sheet)
    undetermined_pairs = pair_undetermined_files(parsed)
    if not undetermined_pairs:
        raise ValueError("no Undetermined R1/R2 FASTQ was found")

    if args.i5_orientation == "auto":
        orientation, orientation_qc = choose_i5_orientation(
            parsed,
            records,
            args.orientation_scan_reads,
            args.max_total_mismatches,
        )
    else:
        orientation = args.i5_orientation
        orientation_qc = {"selected_by_user": 1}
    oriented = orient_expected(records, orientation)
    classifier = FastIndexClassifier(
        oriented,
        args.max_total_mismatches,
        args.max_per_index_mismatches,
        args.minimum_margin,
    )
    selected_compression_backend = args.compression_backend
    if selected_compression_backend == "auto":
        selected_compression_backend = "pigz" if shutil.which("pigz") else "parallel_python"

    detected_assigned_samples = sorted(
        set(item.sample for item in parsed if not item.sample.lower().startswith("undetermined"))
    )
    expected_samples = sorted(record.sample for record in records)
    missing_expected_samples = sorted(set(expected_samples) - set(detected_assigned_samples))
    unexpected_assigned_samples = sorted(set(detected_assigned_samples) - set(expected_samples))
    minimum_pair_distance = minimum_expected_pair_distance(oriented)
    if minimum_pair_distance <= 2 * args.max_total_mismatches:
        warnings.append(
            "Some expected dual-index pairs are close relative to the rescue radius; "
            "review ambiguous counts and do not relax the unique-nearest rule."
        )

    inspection = {
        "input_dir": str(args.input_dir.resolve()),
        "sample_sheet": str(args.sample_sheet.resolve()),
        "detected_fastq_files": len(parsed),
        "detected_samples": sorted(set(item.sample for item in parsed)),
        "expected_samples": expected_samples,
        "missing_expected_samples_in_assigned_fastqs": missing_expected_samples,
        "unexpected_assigned_fastq_samples": unexpected_assigned_samples,
        "undetermined_chunks": len(undetermined_pairs),
        "selected_i5_orientation": orientation,
        "orientation_qc": orientation_qc,
        "minimum_expected_dual_index_distance": minimum_pair_distance,
        "performance": {
            "index_lookup_entries": len(classifier.lookup),
            "requested_compression_backend": args.compression_backend,
            "selected_compression_backend": selected_compression_backend,
            "pigz_threads_per_output": (
                args.pigz_threads_per_file if selected_compression_backend == "pigz" else 0
            ),
        },
        "warnings": warnings,
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "rescue_inspection.json").write_text(
        json.dumps(inspection, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(inspection, indent=2, ensure_ascii=False))
    if args.mode == "inspect":
        return 0

    if missing_expected_samples or unexpected_assigned_samples:
        raise ValueError(
            "SampleSheet Sample_ID values must match assigned FASTQ sample names before rescue. "
            f"Missing={missing_expected_samples}; unexpected={unexpected_assigned_samples}"
        )

    output_fastq_dir = args.outdir / "fastq"
    linked = link_assigned_fastqs(parsed, output_fastq_dir)
    status_counts: Counter[str] = Counter()
    sample_counts: Counter[str] = Counter()
    distance_pairs: Counter[tuple[int, int]] = Counter()
    observed_barcodes: Counter[str] = Counter()
    total_r1_bytes = sum(r1_file.path.stat().st_size for r1_file, _ in undetermined_pairs)
    completed_r1_bytes = 0
    progress = RescueProgressReporter(
        args.outdir / "rescue_progress.json",
        total_r1_bytes,
        args.progress_interval_seconds,
    )
    progress.report(0, 0, 0, 1, len(undetermined_pairs), force=True)
    print(
        f"[rescue] fast index lookup: {len(classifier.lookup):,} entries; "
        f"output compression: {selected_compression_backend}"
        + (
            f" ({args.pigz_threads_per_file} threads/output)"
            if selected_compression_backend == "pigz"
            else ""
        ),
        file=sys.stderr,
        flush=True,
    )

    for chunk_number, (r1_file, r2_file) in enumerate(undetermined_pairs, start=1):
        r1_size = r1_file.path.stat().st_size
        writers = LazyFastqWriters(
            output_fastq_dir,
            chunk_number,
            r2_file is not None,
            args.compresslevel,
            args.compression_backend,
            args.pigz_threads_per_file,
        )
        current_r1_position = [0]
        iterator1 = iter_fastq(
            r1_file.path,
            lambda position: current_r1_position.__setitem__(0, position),
        )
        iterator2 = iter_fastq(r2_file.path) if r2_file else None
        iterator = (
            zip_longest(iterator1, iterator2)
            if iterator2 is not None
            else ((record, None) for record in iterator1)
        )
        try:
            for record1, record2 in iterator:
                if record1 is None or (r2_file is not None and record2 is None):
                    raise ValueError(f"R1/R2 record count mismatch: {r1_file.path}, {r2_file.path}")
                status_counts["total_undetermined_read_pairs_or_reads"] += 1
                processed = status_counts["total_undetermined_read_pairs_or_reads"]
                if processed % args.progress_check_reads == 0:
                    progress.report(
                        processed,
                        sum(sample_counts.values()),
                        completed_r1_bytes + min(current_r1_position[0], r1_size),
                        chunk_number,
                        len(undetermined_pairs),
                    )
                if record2 is not None and canonical_read_id(record1.header) != canonical_read_id(record2.header):
                    status_counts["pair_header_mismatch"] += 1
                    writers.write("Undetermined_residual", record1, record2)
                    continue
                observed1 = extract_header_indexes(record1.header)
                observed2 = extract_header_indexes(record2.header) if record2 is not None else observed1
                if observed1 is not None:
                    observed_barcodes["+".join(observed1)] += 1
                if observed2 is not None and observed2 != observed1:
                    status_counts["paired_header_index_disagreement"] += 1
                    writers.write("Undetermined_residual", record1, record2)
                    continue
                sample, status, _, d7, d5 = classifier.classify(observed1)
                status_counts[status] += 1
                if sample is None:
                    writers.write("Undetermined_residual", record1, record2)
                else:
                    sample_counts[sample] += 1
                    if d7 is not None and d5 is not None:
                        distance_pairs[(d7, d5)] += 1
                    writers.write(sample, record1, record2)
        finally:
            writers.close()
        completed_r1_bytes += r1_size
        progress.report(
            status_counts["total_undetermined_read_pairs_or_reads"],
            sum(sample_counts.values()),
            completed_r1_bytes,
            chunk_number,
            len(undetermined_pairs),
            force=True,
        )

    rescued_total = sum(sample_counts.values())
    total = status_counts["total_undetermined_read_pairs_or_reads"]
    summary_rows: list[dict[str, object]] = []
    for sample in [record.sample for record in oriented]:
        count = sample_counts[sample]
        summary_rows.append(
            {
                "category": "sample",
                "name": sample,
                "read_pairs_or_reads": count,
                "percent_of_undetermined": 100 * count / total if total else 0,
            }
        )
    for status, count in sorted(status_counts.items()):
        summary_rows.append(
            {
                "category": "status",
                "name": status,
                "read_pairs_or_reads": count,
                "percent_of_undetermined": 100 * count / total if total else 0,
            }
        )
    with (args.outdir / "rescue_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    with (args.outdir / "rescue_mismatch_distribution.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=["i7_mismatches", "i5_mismatches", "count"])
        writer.writeheader()
        for (d7, d5), count in sorted(distance_pairs.items()):
            writer.writerow({"i7_mismatches": d7, "i5_mismatches": d5, "count": count})
    with (args.outdir / "top_undetermined_barcodes.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=["barcode", "count"])
        writer.writeheader()
        for barcode, count in observed_barcodes.most_common(1000):
            writer.writerow({"barcode": barcode, "count": count})

    manifest = {
        **inspection,
        "assigned_fastqs_linked": linked,
        "total_undetermined_read_pairs_or_reads": total,
        "rescued_read_pairs_or_reads": rescued_total,
        "rescued_percent": 100 * rescued_total / total if total else 0,
        "sample_rescued_counts": dict(sample_counts),
        "status_counts": dict(status_counts),
        "performance": {
            **inspection["performance"],
            "index_classifier": classifier.stats(),
        },
        "settings": {
            "max_total_mismatches": args.max_total_mismatches,
            "max_per_index_mismatches": args.max_per_index_mismatches,
            "minimum_margin": args.minimum_margin,
            "compresslevel": args.compresslevel,
            "compression_backend": selected_compression_backend,
            "pigz_threads_per_file": args.pigz_threads_per_file,
        },
    }
    (args.outdir / "rescue_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    progress.report(
        total,
        rescued_total,
        total_r1_bytes,
        len(undetermined_pairs),
        len(undetermined_pairs),
        force=True,
        status="completed",
    )
    print(
        f"Rescued {rescued_total:,}/{total:,} Undetermined read pairs/reads "
        f"({manifest['rescued_percent']:.2f}%). Output: {output_fastq_dir}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
