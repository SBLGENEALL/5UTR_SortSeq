#!/usr/bin/env python3
"""Convert NGS_LibraryQC's count matrix into seven-sample Sort-seq inputs."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import pandas as pd


METADATA_COLUMNS = {"sequence_key", "variant_ids", "target_sequence", "length", "gc_percent"}


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t" if path.suffix.lower() in {".tsv", ".txt"} else ",")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--sample-map", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    matrix = read_table(args.matrix)
    mapping = read_table(args.sample_map)
    required_map = {
        "sample_id",
        "ngs_column",
        "sample_type",
        "bin_number",
        "population_fraction",
    }
    missing = required_map.difference(mapping.columns)
    if missing:
        raise ValueError(f"sample map missing columns: {sorted(missing)}")
    if len(mapping) != 7 or set(mapping["sample_type"]) != {"unsorted", "bin"}:
        raise ValueError("sample map must contain one unsorted and six bin rows")
    if mapping["sample_id"].duplicated().any() or mapping["ngs_column"].duplicated().any():
        raise ValueError("sample_id and ngs_column must each be unique")
    missing_matrix = sorted(set(mapping["ngs_column"].astype(str)) - set(matrix.columns))
    if missing_matrix:
        available = [column for column in matrix.columns if column not in METADATA_COLUMNS]
        raise ValueError(
            f"sample-map ngs_column values not found: {missing_matrix}; available sample columns: {available}"
        )
    if "sequence_key" not in matrix.columns and "variant_ids" not in matrix.columns:
        raise ValueError("NGS_LibraryQC matrix must contain sequence_key or variant_ids")
    variant_key = "sequence_key" if "sequence_key" in matrix.columns else "variant_ids"
    variants = pd.DataFrame(
        {
            "variant_id": matrix[variant_key].astype(str),
            "original_variant_ids": matrix.get("variant_ids", matrix[variant_key]).astype(str),
        }
    )
    if "target_sequence" in matrix.columns:
        variants["sequence"] = matrix["target_sequence"].astype(str)
    for column in ["length", "gc_percent"]:
        if column in matrix.columns:
            variants[column] = matrix[column]
    if variants["variant_id"].duplicated().any():
        raise ValueError("variant identifiers in the count matrix are not unique")

    args.outdir.mkdir(parents=True, exist_ok=True)
    counts_dir = args.outdir / "counts"
    counts_dir.mkdir(parents=True, exist_ok=True)
    variants.to_csv(args.outdir / "variants.tsv", sep="\t", index=False)

    sample_rows: list[dict[str, object]] = []
    optional = ["cells_collected", "representative_log10_mfi"]
    for column in optional:
        if column not in mapping.columns:
            mapping[column] = ""
    for _, row in mapping.iterrows():
        sample_id = str(row["sample_id"])
        ngs_column = str(row["ngs_column"])
        count_filename = re.sub(r"[^A-Za-z0-9_.-]+", "_", sample_id) + ".tsv"
        count_table = pd.DataFrame(
            {
                "variant_id": variants["variant_id"],
                "count": pd.to_numeric(matrix[ngs_column], errors="raise").astype(int),
            }
        )
        count_table.to_csv(counts_dir / count_filename, sep="\t", index=False)
        sample_rows.append(
            {
                "sample_id": sample_id,
                "sample_type": row["sample_type"],
                "count_file": f"counts/{count_filename}",
                "id_column": "variant_id",
                "count_column": "count",
                "bin_number": row["bin_number"],
                "population_fraction": row["population_fraction"],
                "cells_collected": row["cells_collected"],
                "representative_log10_mfi": row["representative_log10_mfi"],
            }
        )
    pd.DataFrame(sample_rows).to_csv(args.outdir / "samples.tsv", sep="\t", index=False)

    duplicate_sequences = int(variants["original_variant_ids"].str.contains("|", regex=False).sum())
    with (args.outdir / "conversion_qc.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        writer.writerow(["variants", len(variants)])
        writer.writerow(["collapsed_duplicate_sequence_rows", duplicate_sequences])
        writer.writerow(["matrix_sample_columns", len(matrix.columns) - len(METADATA_COLUMNS.intersection(matrix.columns))])
    print(f"Wrote Sort-seq inputs to {args.outdir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
