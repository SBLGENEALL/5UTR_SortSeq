#!/usr/bin/env python3
"""Create safe configuration files from demultiplexed FASTQs and a run SampleSheet.

The script inspects filenames and reads index metadata. It never moves or
rewrites FASTQ data.
It recognizes names such as ``UTR_bin1_A2_s1_R1_001.fastq.gz`` as sample
``UTR_bin1_A2`` and maps filenames containing bin1..bin6 or unsorted to the
canonical Sort-seq sample IDs.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from rescue_dual_index import ExpectedIndex, discover_fastqs, read_sample_sheet


BIN_FRACTIONS = {1: 0.05, 2: 0.10, 3: 0.15, 4: 0.20, 5: 0.30, 6: 0.20}


def infer_sample_id(fastq_sample: str) -> str | None:
    lower = fastq_sample.lower()
    if lower.startswith("undetermined"):
        return None
    if "unsorted" in lower:
        return "unsorted"
    matches = re.findall(r"(?:^|[^a-z0-9])bin[_-]?([1-6])(?:$|[^0-9])", lower)
    if len(set(matches)) == 1:
        return f"bin{matches[0]}"
    return None


def refuse_existing(paths: list[Path]) -> None:
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing configuration file(s): " + ", ".join(existing)
        )


def normalize_sample_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def find_run_sample_sheet(raw_dir: Path, explicit: Path | None) -> Path | None:
    if explicit is not None:
        if not explicit.is_file():
            raise FileNotFoundError(f"Run SampleSheet not found: {explicit}")
        return explicit.resolve()

    analysis_dir = raw_dir / "Analysis"
    if not analysis_dir.is_dir():
        return None
    candidates = sorted(
        path.resolve()
        for path in analysis_dir.rglob("*.csv")
        if "sample_sheet" in path.name.lower() or "samplesheet" in path.name.lower()
    )
    if len(candidates) > 1:
        raise ValueError(
            "Multiple run SampleSheets were found; select one with --run-sample-sheet: "
            + ", ".join(str(path) for path in candidates)
        )
    return candidates[0] if candidates else None


def match_indexes_to_fastqs(
    run_records: list[ExpectedIndex], assigned_samples: list[str]
) -> dict[str, ExpectedIndex]:
    by_exact = {record.sample: record for record in run_records}
    by_normalized: dict[str, list[ExpectedIndex]] = {}
    for record in run_records:
        by_normalized.setdefault(normalize_sample_name(record.sample), []).append(record)

    matched: dict[str, ExpectedIndex] = {}
    used_source_samples: set[str] = set()
    problems: list[str] = []
    for fastq_sample in assigned_samples:
        record = by_exact.get(fastq_sample)
        if record is None:
            candidates = by_normalized.get(normalize_sample_name(fastq_sample), [])
            if len(candidates) == 1:
                record = candidates[0]
            elif len(candidates) > 1:
                problems.append(
                    f"ambiguous normalized Sample_ID for FASTQ {fastq_sample}: "
                    + ", ".join(candidate.sample for candidate in candidates)
                )
            else:
                problems.append(f"no Sample_ID matches FASTQ prefix {fastq_sample}")
        if record is not None:
            matched[fastq_sample] = record
            used_source_samples.add(record.sample)

    unused = sorted(set(record.sample for record in run_records) - used_source_samples)
    if problems or len(matched) != 7 or unused:
        raise ValueError(
            "The run SampleSheet and assigned FASTQ prefixes must match one-to-one. "
            f"Problems={problems}; unused_Sample_ID={unused}; "
            f"FASTQ_prefixes={assigned_samples}"
        )
    return matched


def write_sample_sheet(
    path: Path,
    ordered_samples: list[str],
    index_by_fastq: dict[str, ExpectedIndex] | None,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        handle.write("[Header]\nFileFormatVersion,2\n\n[Data]\n")
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["Sample_ID", "index", "index2"])
        for sample in ordered_samples:
            if index_by_fastq is None:
                writer.writerow([sample, "REPLACE_I7", "REPLACE_I5"])
            else:
                record = index_by_fastq[sample]
                writer.writerow([sample, record.i7, record.i5])


def write_sample_map(path: Path, canonical_to_fastq: dict[str, str]) -> None:
    fields = [
        "sample_id",
        "ngs_column",
        "sample_type",
        "bin_number",
        "population_fraction",
        "cells_collected",
        "representative_log10_mfi",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for number in range(1, 7):
            canonical = f"bin{number}"
            writer.writerow(
                {
                    "sample_id": canonical,
                    "ngs_column": canonical_to_fastq[canonical],
                    "sample_type": "bin",
                    "bin_number": number,
                    "population_fraction": BIN_FRACTIONS[number],
                    "cells_collected": "",
                    "representative_log10_mfi": "",
                }
            )
        writer.writerow(
            {
                "sample_id": "unsorted",
                "ngs_column": canonical_to_fastq["unsorted"],
                "sample_type": "unsorted",
                "bin_number": "",
                "population_fraction": "",
                "cells_collected": "",
                "representative_log10_mfi": "",
            }
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Detect seven FASTQ sample prefixes and create editable config CSV files"
    )
    parser.add_argument("--raw-dir", required=True, type=Path)
    parser.add_argument("--config-dir", required=True, type=Path)
    parser.add_argument(
        "--run-sample-sheet",
        type=Path,
        help=(
            "Original Illumina/BCL Convert SampleSheet containing Sample_ID, index, "
            "and index2. If omitted, exactly one *sample_sheet*.csv under "
            "RAW_DIR/Analysis is used automatically."
        ),
    )
    args = parser.parse_args()

    parsed, warnings = discover_fastqs(args.raw_dir)
    assigned_samples = sorted(
        {item.sample for item in parsed if not item.sample.lower().startswith("undetermined")}
    )
    if not assigned_samples:
        raise ValueError(f"No assigned R1/R2 FASTQ files were detected under {args.raw_dir}")

    canonical_to_fastq: dict[str, str] = {}
    unclassified: list[str] = []
    duplicates: list[str] = []
    for sample in assigned_samples:
        canonical = infer_sample_id(sample)
        if canonical is None:
            unclassified.append(sample)
        elif canonical in canonical_to_fastq:
            duplicates.append(f"{canonical}: {canonical_to_fastq[canonical]}, {sample}")
        else:
            canonical_to_fastq[canonical] = sample

    expected = [f"bin{x}" for x in range(1, 7)] + ["unsorted"]
    missing = [sample for sample in expected if sample not in canonical_to_fastq]
    if len(assigned_samples) != 7 or unclassified or duplicates or missing:
        raise ValueError(
            "Could not infer exactly one FASTQ prefix for bin1-bin6 and unsorted. "
            f"Detected={assigned_samples}; missing={missing}; unclassified={unclassified}; "
            f"duplicates={duplicates}; filename_warnings={warnings}"
        )

    run_sample_sheet = find_run_sample_sheet(args.raw_dir, args.run_sample_sheet)
    index_by_fastq: dict[str, ExpectedIndex] | None = None
    if run_sample_sheet is not None:
        run_records = read_sample_sheet(run_sample_sheet)
        index_by_fastq = match_indexes_to_fastqs(run_records, assigned_samples)

    args.config_dir.mkdir(parents=True, exist_ok=True)
    detected_path = args.config_dir / "detected_fastq_samples.csv"
    sample_sheet_path = args.config_dir / "SampleSheet.csv"
    sample_map_path = args.config_dir / "sample_map.csv"
    refuse_existing([detected_path, sample_sheet_path, sample_map_path])

    with detected_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["canonical_sample_id", "fastq_sample_prefix"])
        for canonical in expected:
            writer.writerow([canonical, canonical_to_fastq[canonical]])

    ordered_fastq_samples = [canonical_to_fastq[canonical] for canonical in expected]
    write_sample_sheet(sample_sheet_path, ordered_fastq_samples, index_by_fastq)
    write_sample_map(sample_map_path, canonical_to_fastq)

    print(f"Detected sample names: {detected_path}")
    if run_sample_sheet is None:
        print(f"WARNING: no run SampleSheet was found; EDIT i7/i5 before preflight: {sample_sheet_path}")
    else:
        print(f"Imported and validated i7/i5 from: {run_sample_sheet}")
        print(f"Generated rescue SampleSheet: {sample_sheet_path}")
    print(f"Review bin mapping and sorter values: {sample_map_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
