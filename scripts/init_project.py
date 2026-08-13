#!/usr/bin/env python3
"""Create safe starter configuration files from demultiplexed FASTQ names.

The script only inspects filenames. It never reads, moves, or rewrites FASTQ data.
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

from rescue_dual_index import discover_fastqs


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


def write_sample_sheet(path: Path, ordered_samples: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        handle.write("[Header]\nFileFormatVersion,2\n\n[Data]\n")
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["Sample_ID", "index", "index2"])
        for sample in ordered_samples:
            writer.writerow([sample, "REPLACE_I7", "REPLACE_I5"])


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
    write_sample_sheet(sample_sheet_path, ordered_fastq_samples)
    write_sample_map(sample_map_path, canonical_to_fastq)

    print(f"Detected sample names: {detected_path}")
    print(f"EDIT i7/i5 indexes before preflight: {sample_sheet_path}")
    print(f"Review bin mapping and sorter values: {sample_map_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
