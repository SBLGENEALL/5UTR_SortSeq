#!/usr/bin/env python3
"""Export every Sort-seq calculation from a LibraryQC variant-count matrix.

The eight canonical exported calculations are:
  1. raw variant counts
  2. sequencing-depth-normalized within-sample frequencies
  3. FACS-bin-size-corrected population mass
  4. corrected-mass total for each UTR
  5. within-UTR-normalized bin probabilities
  6. conditional High15 probability (primary)
  7. per-bin score contributions and expected bin score (supporting)
  8. relative-enrichment profile p_ib / w_b (visualization)

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
        description="Export eight auditable Sort-seq calculation steps as CSV files."
    )
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--sample-map", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--pseudocount", type=float, default=0.5)
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
    fractions = pd.to_numeric(
        bins["population_fraction"], errors="raise"
    ).to_numpy(float)
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
    if args.pseudocount <= 0:
        raise ValueError("pseudocount must be positive")
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
    n_variants = len(metadata)
    smooth_bins = (raw_bins + args.pseudocount).div(
        depths[bin_ids] + args.pseudocount * n_variants,
        axis=1,
    )
    smooth_unsorted = (raw_unsorted + args.pseudocount) / (
        depths[unsorted_id] + args.pseudocount * n_variants
    )
    top15_combined_frequency = (
        smooth_bins.iloc[:, :2].mul(fractions.iloc[:2], axis=1).sum(axis=1)
        / high15_fraction
    )
    bin1_vs_unsorted_enrichment = smooth_bins.iloc[:, 0] / smooth_unsorted
    bin2_vs_unsorted_enrichment = smooth_bins.iloc[:, 1] / smooth_unsorted
    top15_vs_unsorted_enrichment = top15_combined_frequency / smooth_unsorted
    step5 = metadata.copy()
    for number, sample_id in enumerate(bin_ids, start=1):
        step5[f"bin{number}_score_weight"] = float(weights[sample_id])
        step5[f"bin{number}_score_contribution"] = contributions[sample_id].to_numpy()
    step5["expected_bin_score"] = expected_score.to_numpy()
    step5["neutral_baseline_score"] = neutral_baseline
    step5["score_shift_from_neutral"] = expected_score.to_numpy() - neutral_baseline
    step5["high15_probability"] = high15_probability.to_numpy()
    step5["high15_enrichment"] = high15_probability.to_numpy() / high15_fraction
    step5["bin1_vs_unsorted_enrichment"] = bin1_vs_unsorted_enrichment.to_numpy()
    step5["bin2_vs_unsorted_enrichment"] = bin2_vs_unsorted_enrichment.to_numpy()
    step5["bin1_vs_unsorted_log2_enrichment"] = np.log2(
        bin1_vs_unsorted_enrichment.to_numpy()
    )
    step5["bin2_vs_unsorted_log2_enrichment"] = np.log2(
        bin2_vs_unsorted_enrichment.to_numpy()
    )
    step5["top15_combined_frequency"] = top15_combined_frequency.to_numpy()
    step5["top15_vs_unsorted_enrichment"] = top15_vs_unsorted_enrichment.to_numpy()
    step5["top15_vs_unsorted_log2_enrichment"] = np.log2(
        top15_vs_unsorted_enrichment.to_numpy()
    )

    # Canonical Step 6: primary conditional high-tail endpoint.
    step6 = metadata.copy()
    step6["step6_high15_probability"] = high15_probability.to_numpy()
    step6["high15_probability"] = high15_probability.to_numpy()
    step6["high15_nominal_population_fraction"] = high15_fraction
    step6["step8_high15_relative_enrichment"] = (
        high15_probability.to_numpy() / high15_fraction
    )
    step6["step8_high15_log2_relative_enrichment"] = np.log2(
        np.clip(step6["step8_high15_relative_enrichment"].to_numpy(), 1e-12, None)
    )

    # Canonical Step 7: supporting whole-distribution ordinal score.
    step7 = metadata.copy()
    for number, sample_id in enumerate(bin_ids, start=1):
        step7[f"bin{number}_score_weight"] = float(weights[sample_id])
        step7[f"bin{number}_score_contribution"] = contributions[
            sample_id
        ].to_numpy()
    step7["step7_expected_bin_score"] = expected_score.to_numpy()
    step7["expected_bin_score"] = expected_score.to_numpy()
    step7["neutral_baseline_score"] = neutral_baseline
    step7["score_shift_from_neutral"] = expected_score.to_numpy() - neutral_baseline

    # Canonical Step 8: p_ib / w_b. A neutral UTR equals one in every bin;
    # log2 values equal zero. This is for profile shape, not a second probability.
    relative_enrichment = probability.div(fractions, axis=1)
    step8 = metadata.copy()
    for number, sample_id in enumerate(bin_ids, start=1):
        value = relative_enrichment[sample_id]
        step8[f"bin{number}_relative_enrichment"] = value.to_numpy()
        step8[f"bin{number}_log2_relative_enrichment"] = np.log2(
            np.clip(value.to_numpy(), 1e-12, None)
        )
    step8["step8_high15_relative_enrichment"] = (
        high15_probability.to_numpy() / high15_fraction
    )
    step8["step8_high15_log2_relative_enrichment"] = np.log2(
        np.clip(step8["step8_high15_relative_enrichment"].to_numpy(), 1e-12, None)
    )

    reference_id = resolve_reference_id(metadata, args.reference_variant_id)
    if reference_id is not None:
        score_by_id = pd.Series(expected_score.to_numpy(), index=metadata["variant_id"])
        high15_by_id = pd.Series(high15_probability.to_numpy(), index=metadata["variant_id"])
        reference_score = float(score_by_id.loc[reference_id])
        reference_high15 = float(high15_by_id.loc[reference_id])
        top15_unsorted_by_id = pd.Series(
            top15_vs_unsorted_enrichment.to_numpy(), index=metadata["variant_id"]
        )
        reference_top15_unsorted = float(top15_unsorted_by_id.loc[reference_id])
        step5["is_reference_variant"] = metadata["variant_id"].eq(reference_id)
        step5["reference_expected_bin_score"] = reference_score
        step5["delta_score_vs_reference"] = expected_score.to_numpy() - reference_score
        step5["reference_high15_probability"] = reference_high15
        step5["delta_high15_probability_vs_reference"] = (
            high15_probability.to_numpy() - reference_high15
        )
        step5["high15_fold_vs_reference"] = (
            high15_probability.to_numpy() / reference_high15
        )
        step5["high15_log2_fold_vs_reference"] = np.log2(
            high15_probability.to_numpy() / reference_high15
        )
        step5["reference_top15_vs_unsorted_enrichment"] = reference_top15_unsorted
        step5["delta_top15_log2_enrichment_vs_reference"] = np.log2(
            top15_vs_unsorted_enrichment.to_numpy() / reference_top15_unsorted
        )
        step5["top15_enrichment_fold_vs_reference"] = (
            top15_vs_unsorted_enrichment.to_numpy() / reference_top15_unsorted
        )

    # A single wide audit table makes one-row tracing convenient.
    audit = step1.copy()
    for table in [step2, step3, step4, step5, step6, step7, step8]:
        extra = table.drop(columns=list(metadata.columns), errors="ignore")
        new_columns = [column for column in extra.columns if column not in audit.columns]
        audit = pd.concat(
            [audit.reset_index(drop=True), extra[new_columns].reset_index(drop=True)],
            axis=1,
        )
    audit["min200_exploratory_pass"] = (
        (audit["unsorted_raw_count"] >= 50)
        & (audit["total_6bin_raw_count"] >= 200)
        & (audit["detected_in_n_bins"] >= 3)
    )
    audit["maximum_bin_probability"] = audit[
        [f"bin{x}_probability" for x in range(1, 7)]
    ].max(axis=1)
    audit["single_bin_dominance_flag"] = audit["maximum_bin_probability"] >= 0.85
    audit["top15_read_support_pass"] = (
        (audit["total_6bin_raw_count"] >= 200)
        & (audit["bin1_plus_bin2_raw_count"] >= 20)
    )
    audit["top15_unsorted_qc_pass"] = audit["unsorted_raw_count"] >= 50
    audit["top15_both_bins_enriched"] = (
        (audit["bin1_vs_unsorted_enrichment"] >= 1.0)
        & (audit["bin2_vs_unsorted_enrichment"] >= 1.0)
    )
    audit["high15_primary_rank"] = audit["high15_probability"].where(
        audit["top15_read_support_pass"]
    ).rank(ascending=False, method="average")
    audit["top15_rank"] = audit["high15_primary_rank"]
    audit = audit.sort_values("expected_bin_score", ascending=False, na_position="last")

    top15_columns = [
        "variant_id",
        "high15_primary_rank",
        "top15_rank",
        "high15_probability",
        "high15_enrichment",
        "expected_bin_score",
        "top15_vs_unsorted_enrichment",
        "top15_vs_unsorted_log2_enrichment",
        "bin1_vs_unsorted_enrichment",
        "bin2_vs_unsorted_enrichment",
        "unsorted_raw_count",
        "total_6bin_raw_count",
        "bin1_plus_bin2_raw_count",
        "top15_both_bins_enriched",
        "single_bin_dominance_flag",
    ]
    if reference_id is not None:
        top15_columns.extend(
            [
                "delta_top15_log2_enrichment_vs_reference",
                "top15_enrichment_fold_vs_reference",
                "delta_high15_probability_vs_reference",
                "high15_fold_vs_reference",
                "high15_log2_fold_vs_reference",
            ]
        )
    top15_ranking = audit.loc[
        audit["top15_read_support_pass"], top15_columns
    ].sort_values("high15_primary_rank", ascending=True)

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
            {
                "check": "maximum_abs_step8_high15_identity_error",
                "value": float(
                    np.nanmax(
                        np.abs(
                            step8["step8_high15_relative_enrichment"].to_numpy()
                            - high15_probability.to_numpy() / high15_fraction
                        )
                    )
                ),
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
    write_csv(top15_ranking, args.outdir / "08_high15_primary_ranking.csv")
    write_csv(
        top15_ranking.sort_values(
            "top15_vs_unsorted_log2_enrichment", ascending=False
        ),
        args.outdir / "09_top15_unsorted_enrichment_secondary.csv",
    )

    # Canonical, descriptive filenames. Legacy numbered outputs above are kept
    # so previous notebooks and reports continue to work unchanged.
    step4_total = metadata.copy()
    step4_total["sum_6bin_population_mass"] = corrected_mass.sum(axis=1).to_numpy()
    write_csv(step1, args.outdir / "step01_raw_counts.csv")
    write_csv(step2, args.outdir / "step02_depth_normalized_frequency.csv")
    write_csv(step3, args.outdir / "step03_bin_size_corrected_mass.csv")
    write_csv(step4_total, args.outdir / "step04_corrected_mass_total.csv")
    write_csv(step4, args.outdir / "step05_within_utr_probability.csv")
    write_csv(step6, args.outdir / "step06_high15_probability.csv")
    write_csv(step7, args.outdir / "step07_expected_score.csv")
    write_csv(step8, args.outdir / "step08_relative_enrichment_profile.csv")
    write_csv(audit, args.outdir / "step09_all_steps_audit.csv")

    manifest = {
        "matrix": str(args.matrix.resolve()),
        "sample_map": str(args.sample_map.resolve()),
        "variants": len(metadata),
        "population_fraction_source": fraction_source,
        "neutral_baseline_score": neutral_baseline,
        "reference_variant_id": reference_id,
        "formula": "p_ib = w_b*(c_ib/N_b) / sum_k[w_k*(c_ik/N_k)]",
        "score": "S_i = sum_b[p_ib*(7-bin_number_b)]",
        "primary_endpoint": "H_i = p_i1 + p_i2",
        "relative_enrichment_profile": "R_ib = p_ib / w_b; L_ib = log2(R_ib)",
        "step8_scalar_identity": "H_i/(w1+w2), identical ranking to Step 6",
        "primary_rank": "high15_probability descending (bootstrap is added by analyze_sortseq.py)",
        "top15_unsorted_enrichment": (
            "secondary E_i = [(w1*f_i1+w2*f_i2)/(w1+w2)] / f_i_unsorted"
        ),
        "top15_pseudocount": args.pseudocount,
    }
    (args.outdir / "scoring_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"Eight scoring-step CSVs written to: {args.outdir}")
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
