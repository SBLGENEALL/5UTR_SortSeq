#!/usr/bin/env python3
"""Export every Sort-seq scoring step from a LibraryQC variant-count matrix.

The five exported calculations are:
  1. raw variant counts
  2. sequencing-depth-normalized within-sample frequencies
  3. FACS-bin-size-corrected population mass
  4. within-UTR-normalized bin probabilities
  5. per-bin score contributions and final expected bin score

Unsorted is exported for coverage/QC but is not included in the six-bin
fluorescence score.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


METADATA_COLUMNS = {
    "sequence_key",
    "variant_ids",
    "target_sequence",
    "length",
    "gc_percent",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export five auditable Sort-seq scoring steps as CSV files."
    )
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--sample-map", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument(
        "--reference-variant-id",
        default="auto",
        help="Reference ID; auto recognizes original/orginal, and none disables it.",
    )
    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    separator = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=separator)


def write_csv(table: pd.DataFrame, path: Path) -> None:
    # utf-8-sig opens cleanly in Excel while remaining valid UTF-8 CSV.
    table.to_csv(path, index=False, encoding="utf-8-sig")


def resolve_reference_id(metadata: pd.DataFrame, requested: str) -> str | None:
    value = str(requested).strip()
    if value.lower() in {"", "none", "off", "false"}:
        return None
    aliases = {"original", "orginal"} if value.lower() == "auto" else {value.lower()}
    direct = metadata["variant_id"].astype(str).str.strip().str.lower().isin(aliases)
    original = metadata["original_variant_ids"].fillna("").astype(str).map(
        lambda item: any(part.strip().lower() in aliases for part in item.split("|"))
    )
    matches = metadata.loc[direct | original, "variant_id"].astype(str).tolist()
    if not matches:
        if value.lower() == "auto":
            return None
        raise ValueError(f"reference variant ID not found: {value!r}")
    if len(matches) > 1:
        raise ValueError(f"reference variant is ambiguous: {matches}")
    return matches[0]


def load_inputs(matrix_path: Path, sample_map_path: Path):
    matrix = read_table(matrix_path)
    mapping = read_table(sample_map_path)
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
    mapping = mapping.copy()
    mapping["sample_id"] = mapping["sample_id"].astype(str)
    mapping["ngs_column"] = mapping["ngs_column"].astype(str)
    mapping["sample_type"] = mapping["sample_type"].astype(str).str.lower()
    if len(mapping) != 7:
        raise ValueError("sample map must have exactly six bins and one unsorted row")
    if len(mapping[mapping["sample_type"].eq("bin")]) != 6:
        raise ValueError("sample map must have exactly six bin rows")
    if len(mapping[mapping["sample_type"].eq("unsorted")]) != 1:
        raise ValueError("sample map must have exactly one unsorted row")
    missing_columns = sorted(set(mapping["ngs_column"]) - set(matrix.columns))
    if missing_columns:
        raise ValueError(f"matrix is missing mapped NGS columns: {missing_columns}")

    variant_key = "sequence_key" if "sequence_key" in matrix.columns else "variant_ids"
    if variant_key not in matrix.columns:
        raise ValueError("matrix must contain sequence_key or variant_ids")
    metadata = pd.DataFrame(
        {
            "variant_id": matrix[variant_key].astype(str),
            "original_variant_ids": matrix.get("variant_ids", matrix[variant_key]).astype(str),
        }
    )
    for source, target in [
        ("target_sequence", "target_sequence"),
        ("length", "length"),
        ("gc_percent", "gc_percent"),
    ]:
        if source in matrix.columns:
            metadata[target] = matrix[source]
    if metadata["variant_id"].duplicated().any():
        raise ValueError("variant IDs in the matrix must be unique")

    counts = pd.DataFrame(index=metadata["variant_id"])
    for _, row in mapping.iterrows():
        counts[str(row["sample_id"])] = pd.to_numeric(
            matrix[str(row["ngs_column"])], errors="raise"
        ).to_numpy(dtype=float)
    if (counts < 0).any().any():
        raise ValueError("counts cannot be negative")

    bins = mapping[mapping["sample_type"].eq("bin")].copy()
    bins["bin_number"] = pd.to_numeric(bins["bin_number"], errors="raise").astype(int)
    bins = bins.sort_values("bin_number")
    if bins["bin_number"].tolist() != [1, 2, 3, 4, 5, 6]:
        raise ValueError("bin_number must contain 1,2,3,4,5,6 exactly once")
    if "cells_collected" in bins.columns:
        collected = pd.to_numeric(bins["cells_collected"], errors="coerce")
    else:
        collected = pd.Series(np.nan, index=bins.index)
    if collected.notna().all() and (collected > 0).all():
        fractions = collected.to_numpy(dtype=float) / collected.sum()
        fraction_source = "cells_collected"
    elif collected.notna().any():
        raise ValueError("cells_collected must be positive for all six bins or blank for all")
    else:
        fractions = pd.to_numeric(bins["population_fraction"], errors="raise").to_numpy(float)
        if fractions.sum() > 1.5:
            fractions = fractions / 100.0
        if (fractions <= 0).any() or not np.isclose(fractions.sum(), 1.0, atol=0.02):
            raise ValueError("six population fractions must be positive and sum to 1 or 100")
        fractions = fractions / fractions.sum()
        fraction_source = "population_fraction"
    bins["effective_population_fraction"] = fractions
    bins["score_weight"] = 7 - bins["bin_number"]
    return matrix, metadata, mapping, bins, counts, fraction_source


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    _, metadata, mapping, bins, counts, fraction_source = load_inputs(
        args.matrix, args.sample_map
    )
    bin_ids = bins["sample_id"].astype(str).tolist()
    unsorted_id = str(
        mapping.loc[mapping["sample_type"].eq("unsorted"), "sample_id"].iloc[0]
    )
    depths = counts.sum(axis=0)
    if (depths <= 0).any():
        raise ValueError(f"sample depth is zero: {depths[depths <= 0].index.tolist()}")

    fractions = pd.Series(
        bins["effective_population_fraction"].to_numpy(float), index=bin_ids
    )
    weights = pd.Series(bins["score_weight"].to_numpy(float), index=bin_ids)
    raw_bins = counts[bin_ids]
    raw_unsorted = counts[unsorted_id]

    # Step 1: raw counts.
    step1 = metadata.copy()
    for number, sample_id in enumerate(bin_ids, start=1):
        step1[f"bin{number}_raw_count"] = raw_bins[sample_id].to_numpy()
    step1["unsorted_raw_count"] = raw_unsorted.to_numpy()
    step1["total_6bin_raw_count"] = raw_bins.sum(axis=1).to_numpy()
    step1["bin1_plus_bin2_raw_count"] = raw_bins.iloc[:, :2].sum(axis=1).to_numpy()
    step1["detected_in_n_bins"] = (raw_bins > 0).sum(axis=1).to_numpy()

    # Step 2: P(UTR | bin), correcting unequal NGS sample depth.
    depth_frequency = raw_bins.div(depths[bin_ids], axis=1)
    unsorted_frequency = raw_unsorted / depths[unsorted_id]
    step2 = metadata.copy()
    for number, sample_id in enumerate(bin_ids, start=1):
        step2[f"bin{number}_depth_normalized_frequency"] = depth_frequency[
            sample_id
        ].to_numpy()
    step2["unsorted_depth_normalized_frequency"] = unsorted_frequency.to_numpy()

    # Step 3: estimated joint population mass P(UTR, bin), up to a common scale.
    corrected_mass = depth_frequency.mul(fractions, axis=1)
    step3 = metadata.copy()
    for number, sample_id in enumerate(bin_ids, start=1):
        step3[f"bin{number}_population_mass"] = corrected_mass[sample_id].to_numpy()
    step3["sum_6bin_population_mass"] = corrected_mass.sum(axis=1).to_numpy()

    # Step 4: P(bin | UTR, target gate). Rows sum to one when mass is nonzero.
    mass_sum = corrected_mass.sum(axis=1).replace(0, np.nan)
    probability = corrected_mass.div(mass_sum, axis=0)
    step4 = metadata.copy()
    for number, sample_id in enumerate(bin_ids, start=1):
        step4[f"bin{number}_probability"] = probability[sample_id].to_numpy()
    step4["probability_sum_check"] = probability.sum(axis=1, min_count=1).to_numpy()

    # Step 5: contribution of each bin to the expected ordinal fluorescence score.
    contributions = probability.mul(weights, axis=1)
    expected_score = contributions.sum(axis=1, min_count=1)
    neutral_baseline = float(np.dot(fractions.to_numpy(), weights.to_numpy()))
    high15_probability = probability.iloc[:, :2].sum(axis=1, min_count=1)
    high15_fraction = float(fractions.iloc[:2].sum())
    step5 = metadata.copy()
    for number, sample_id in enumerate(bin_ids, start=1):
        step5[f"bin{number}_score_weight"] = float(weights[sample_id])
        step5[f"bin{number}_score_contribution"] = contributions[sample_id].to_numpy()
    step5["expected_bin_score"] = expected_score.to_numpy()
    step5["neutral_baseline_score"] = neutral_baseline
    step5["score_shift_from_neutral"] = expected_score.to_numpy() - neutral_baseline
    step5["high15_probability"] = high15_probability.to_numpy()
    step5["high15_enrichment"] = high15_probability.to_numpy() / high15_fraction

    reference_id = resolve_reference_id(metadata, args.reference_variant_id)
    if reference_id is not None:
        score_by_id = pd.Series(expected_score.to_numpy(), index=metadata["variant_id"])
        high15_by_id = pd.Series(high15_probability.to_numpy(), index=metadata["variant_id"])
        reference_score = float(score_by_id.loc[reference_id])
        reference_high15 = float(high15_by_id.loc[reference_id])
        step5["is_reference_variant"] = metadata["variant_id"].eq(reference_id)
        step5["reference_expected_bin_score"] = reference_score
        step5["delta_score_vs_reference"] = expected_score.to_numpy() - reference_score
        step5["reference_high15_probability"] = reference_high15
        step5["delta_high15_probability_vs_reference"] = (
            high15_probability.to_numpy() - reference_high15
        )

    # A single wide audit table makes one-row tracing convenient.
    audit = step1.copy()
    for table in [step2, step3, step4, step5]:
        extra = table.drop(columns=list(metadata.columns), errors="ignore")
        audit = pd.concat([audit.reset_index(drop=True), extra.reset_index(drop=True)], axis=1)
    audit["min200_exploratory_pass"] = (
        (audit["unsorted_raw_count"] >= 50)
        & (audit["total_6bin_raw_count"] >= 200)
        & (audit["detected_in_n_bins"] >= 3)
    )
    audit["maximum_bin_probability"] = audit[
        [f"bin{x}_probability" for x in range(1, 7)]
    ].max(axis=1)
    audit["single_bin_dominance_flag"] = audit["maximum_bin_probability"] >= 0.85
    audit = audit.sort_values("expected_bin_score", ascending=False, na_position="last")

    parameter_rows = []
    mapped_by_id = mapping.set_index("sample_id")
    for _, row in bins.iterrows():
        sample_id = str(row["sample_id"])
        parameter_rows.append(
            {
                "sample_id": sample_id,
                "ngs_column": mapped_by_id.loc[sample_id, "ngs_column"],
                "sample_type": "bin",
                "bin_number": int(row["bin_number"]),
                "total_assigned_utr_reads": float(depths[sample_id]),
                "population_fraction": float(row["effective_population_fraction"]),
                "score_weight": float(row["score_weight"]),
                "population_fraction_source": fraction_source,
            }
        )
    parameter_rows.append(
        {
            "sample_id": unsorted_id,
            "ngs_column": mapped_by_id.loc[unsorted_id, "ngs_column"],
            "sample_type": "unsorted",
            "bin_number": "",
            "total_assigned_utr_reads": float(depths[unsorted_id]),
            "population_fraction": "",
            "score_weight": "",
            "population_fraction_source": "not_used_in_score",
        }
    )
    parameters = pd.DataFrame(parameter_rows)
    checks = pd.DataFrame(
        [
            {"check": "variants", "value": len(metadata)},
            {"check": "neutral_baseline_score", "value": neutral_baseline},
            {"check": "high15_population_fraction", "value": high15_fraction},
            {
                "check": "maximum_abs_probability_sum_error",
                "value": float((step4["probability_sum_check"].dropna() - 1).abs().max()),
            },
            {"check": "reference_variant_id", "value": reference_id or "not_detected"},
        ]
    )

    write_csv(parameters, args.outdir / "00_sample_parameters.csv")
    write_csv(step1, args.outdir / "01_raw_counts.csv")
    write_csv(step2, args.outdir / "02_depth_normalized_frequency.csv")
    write_csv(step3, args.outdir / "03_bin_size_corrected_mass.csv")
    write_csv(step4, args.outdir / "04_within_utr_bin_probability.csv")
    write_csv(step5, args.outdir / "05_score_contributions_and_final_score.csv")
    write_csv(audit, args.outdir / "06_all_steps_combined_audit.csv")
    write_csv(checks, args.outdir / "07_calculation_checks.csv")

    manifest = {
        "matrix": str(args.matrix.resolve()),
        "sample_map": str(args.sample_map.resolve()),
        "variants": len(metadata),
        "population_fraction_source": fraction_source,
        "neutral_baseline_score": neutral_baseline,
        "reference_variant_id": reference_id,
        "formula": "p_ib = w_b*(c_ib/N_b) / sum_k[w_k*(c_ik/N_k)]",
        "score": "S_i = sum_b[p_ib*(7-bin_number_b)]",
    }
    (args.outdir / "scoring_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"Scoring-step CSVs written to: {args.outdir}")
    print(f"Variants: {len(metadata):,}")
    print(f"Neutral baseline score: {neutral_baseline:.6f}")
    print(f"Reference: {reference_id or 'not detected'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
