#!/usr/bin/env python3
"""Compare Sort-seq Steps 6, 7, and 8 on one read-supported UTR set.

Step 6 is the primary conditional High15 probability, Step 7 is the
supporting six-bin ordinal score, and Step 8 is the relative-enrichment
representation p_ib / w_b. The scalar Step-8 High15 enrichment is exactly
Step6 / (w1 + w2), so its rank must be identical to Step 6. This script
checks that identity and treats any disagreement as a calculation error.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, pearsonr, spearmanr


DEFAULT_POPULATION_FRACTIONS = np.asarray(
    [0.05, 0.10, 0.15, 0.20, 0.30, 0.20], dtype=float
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare Step 6 High15 probability, Step 7 expected-bin score, "
            "and Step 8 relative enrichment on the same UTR universe."
        )
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument(
        "--sample-map",
        type=Path,
        help="Optional sample map used to read the six nominal population fractions.",
    )
    parser.add_argument("--min-total-count", type=int, default=200)
    parser.add_argument(
        "--min-unsorted-count",
        type=int,
        default=0,
        help=(
            "Optional unsorted QC cutoff. It is not part of the Step 6/7/8 formulas; "
            "the default of zero leaves the conditional phenotype universe unchanged."
        ),
    )
    parser.add_argument("--min-high-bin-count", type=int, default=20)
    parser.add_argument("--top-n", type=int, default=50)
    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    separator = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=separator)


def read_population_fractions(sample_map: Path | None) -> np.ndarray:
    if sample_map is None:
        return DEFAULT_POPULATION_FRACTIONS.copy()
    mapping = read_table(sample_map)
    required = {"sample_type", "bin_number", "population_fraction"}
    missing = sorted(required.difference(mapping.columns))
    if missing:
        raise ValueError(f"sample map is missing columns: {missing}")
    bins = mapping[
        mapping["sample_type"].astype(str).str.lower().eq("bin")
    ].copy()
    bins["bin_number"] = pd.to_numeric(bins["bin_number"], errors="raise").astype(int)
    bins = bins.sort_values("bin_number")
    if bins["bin_number"].tolist() != [1, 2, 3, 4, 5, 6]:
        raise ValueError("sample map must contain bin_number 1 through 6 exactly once")
    fractions = pd.to_numeric(
        bins["population_fraction"], errors="raise"
    ).to_numpy(dtype=float)
    if fractions.sum() > 1.5:
        fractions = fractions / 100.0
    if (fractions <= 0).any() or not np.isclose(fractions.sum(), 1.0, atol=0.02):
        raise ValueError("population fractions must be positive and sum to 1 or 100")
    return fractions / fractions.sum()


def as_bool(values: pd.Series) -> pd.Series:
    return values.astype(str).str.lower().isin({"true", "t", "1"})


def safe_correlation(
    x: pd.Series, y: pd.Series, method: str
) -> tuple[float, float, int]:
    finite = np.isfinite(x) & np.isfinite(y)
    left = x[finite]
    right = y[finite]
    n = int(finite.sum())
    if n < 3 or left.nunique() < 2 or right.nunique() < 2:
        return float("nan"), float("nan"), n
    if method == "spearman":
        result = spearmanr(left, right)
    elif method == "pearson":
        result = pearsonr(left, right)
    elif method == "kendall":
        result = kendalltau(left, right)
    else:
        raise ValueError(f"unknown correlation method: {method}")
    return float(result.statistic), float(result.pvalue), n


def rank_overlap(
    table: pd.DataFrame,
    top_n: int,
    left_rank: str = "score_rank_filtered",
    right_rank: str = "top15_rank_filtered",
    pair: str = "step6_vs_step7",
) -> dict[str, object]:
    actual_n = min(top_n, len(table))
    if actual_n == 0:
        return {
            "pair": pair,
            "requested_top_n": top_n,
            "actual_top_n": 0,
            "overlap_count": 0,
            "overlap_percent_of_each_list": float("nan"),
            "jaccard": float("nan"),
        }
    left_ids = set(table.nsmallest(actual_n, left_rank)["variant_id"])
    right_ids = set(table.nsmallest(actual_n, right_rank)["variant_id"])
    overlap = left_ids & right_ids
    union = left_ids | right_ids
    return {
        "pair": pair,
        "requested_top_n": top_n,
        "actual_top_n": actual_n,
        "overlap_count": len(overlap),
        "overlap_percent_of_each_list": 100 * len(overlap) / actual_n,
        "jaccard": len(overlap) / len(union) if union else float("nan"),
    }


def analyze_metric_comparison(
    result: pd.DataFrame,
    min_total_count: int = 200,
    min_unsorted_count: int = 0,
    min_high_bin_count: int = 20,
    top_n: int = 50,
    population_fractions: np.ndarray | list[float] | None = None,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    dict[str, object],
]:
    fractions = np.asarray(
        DEFAULT_POPULATION_FRACTIONS
        if population_fractions is None
        else population_fractions,
        dtype=float,
    )
    if fractions.shape != (6,) or (fractions <= 0).any():
        raise ValueError("population_fractions must contain six positive values")
    fractions = fractions / fractions.sum()
    high15_fraction = float(fractions[:2].sum())

    probability_columns = [f"bin{x}_probability" for x in range(1, 7)]
    required = {
        "variant_id",
        "total_6bin_count",
        "unsorted_count",
        "high_bin_raw_count",
        "expected_bin_score",
        "high15_probability",
        *probability_columns,
    }
    missing = sorted(required.difference(result.columns))
    if missing:
        raise ValueError(f"input is missing columns: {missing}")
    if min_total_count < 0 or min_unsorted_count < 0 or min_high_bin_count < 0:
        raise ValueError("read-count cutoffs cannot be negative")
    if top_n <= 0:
        raise ValueError("top_n must be positive")

    output = result.copy()
    numeric_columns = [
        "total_6bin_count",
        "unsorted_count",
        "high_bin_raw_count",
        "expected_bin_score",
        "high15_probability",
        *probability_columns,
    ]
    if "top15_vs_unsorted_log2_enrichment" in output.columns:
        numeric_columns.append("top15_vs_unsorted_log2_enrichment")
    for column in numeric_columns:
        output[column] = pd.to_numeric(output[column], errors="coerce")

    # Canonical Step 6/7/8 columns. Step 8 is both a six-bin profile and a
    # scalar high-tail summary. The scalar is a fixed transform of Step 6.
    output["step6_high15_probability"] = output["high15_probability"]
    output["step7_expected_bin_score"] = output["expected_bin_score"]
    for index, fraction in enumerate(fractions, start=1):
        relative = output[f"bin{index}_probability"] / float(fraction)
        output[f"bin{index}_relative_enrichment"] = relative
        output[f"bin{index}_log2_relative_enrichment"] = np.log2(
            relative.clip(lower=1e-12)
        )
    output["step8_high15_relative_enrichment"] = (
        output["step6_high15_probability"] / high15_fraction
    )
    output["step8_high15_log2_relative_enrichment"] = np.log2(
        output["step8_high15_relative_enrichment"].clip(lower=1e-12)
    )

    finite_metrics = (
        np.isfinite(output["step6_high15_probability"])
        & np.isfinite(output["step7_expected_bin_score"])
        & np.isfinite(output["step8_high15_relative_enrichment"])
    )
    output["total_count_filter_pass"] = (
        output["total_6bin_count"] >= min_total_count
    )
    output["unsorted_qc_filter_pass"] = (
        output["unsorted_count"] >= min_unsorted_count
    )
    output["high_bin_support_filter_pass"] = (
        output["high_bin_raw_count"] >= min_high_bin_count
    )
    output["comparison_read_support_pass"] = (
        output["total_count_filter_pass"]
        & output["unsorted_qc_filter_pass"]
        & output["high_bin_support_filter_pass"]
        & finite_metrics
    )

    supported = output.loc[output["comparison_read_support_pass"]].copy()
    supported["step6_rank_filtered"] = supported[
        "step6_high15_probability"
    ].rank(method="min", ascending=False)
    supported["step7_rank_filtered"] = supported[
        "step7_expected_bin_score"
    ].rank(method="min", ascending=False)
    supported["step8_rank_filtered"] = supported[
        "step8_high15_relative_enrichment"
    ].rank(method="min", ascending=False)

    # Backward-compatible names used by existing reports and R figures.
    supported["top15_rank_filtered"] = supported["step6_rank_filtered"]
    supported["score_rank_filtered"] = supported["step7_rank_filtered"]
    supported["score_rank_percentile"] = supported[
        "step7_expected_bin_score"
    ].rank(method="average", pct=True)
    supported["top15_rank_percentile"] = supported[
        "step6_high15_probability"
    ].rank(method="average", pct=True)
    supported["score_rank_minus_top15_rank"] = (
        supported["step7_rank_filtered"] - supported["step6_rank_filtered"]
    )
    supported["step6_rank_minus_step8_rank"] = (
        supported["step6_rank_filtered"] - supported["step8_rank_filtered"]
    )

    actual_top_n = min(top_n, len(supported))
    supported["top_score_list_flag"] = False
    supported["top15_list_flag"] = False
    supported["top_step8_list_flag"] = False
    if actual_top_n:
        supported.loc[
            supported.nsmallest(actual_top_n, "step7_rank_filtered").index,
            "top_score_list_flag",
        ] = True
        supported.loc[
            supported.nsmallest(actual_top_n, "step6_rank_filtered").index,
            "top15_list_flag",
        ] = True
        supported.loc[
            supported.nsmallest(actual_top_n, "step8_rank_filtered").index,
            "top_step8_list_flag",
        ] = True
    supported["top_list_membership"] = np.select(
        [
            supported["top_score_list_flag"] & supported["top15_list_flag"],
            supported["top_score_list_flag"],
            supported["top15_list_flag"],
        ],
        ["consensus", "score_only", "top15_only"],
        default="neither",
    )

    if "is_reference_variant" in supported.columns:
        is_reference = as_bool(supported["is_reference_variant"])
    else:
        is_reference = supported["variant_id"].astype(str).str.lower().isin(
            {"original", "orginal"}
        )
    supported["is_reference_variant"] = is_reference
    reference = supported.loc[is_reference]
    if len(reference) > 1:
        raise ValueError("more than one reference variant was identified")
    if len(reference) == 1:
        reference_score = float(reference["step7_expected_bin_score"].iloc[0])
        reference_top15 = float(reference["step6_high15_probability"].iloc[0])
        supported["score_above_reference_filtered"] = (
            supported["step7_expected_bin_score"] > reference_score
        )
        supported["top15_above_reference_filtered"] = (
            supported["step6_high15_probability"] > reference_top15
        )
        supported["reference_quadrant"] = np.select(
            [
                supported["score_above_reference_filtered"]
                & supported["top15_above_reference_filtered"],
                supported["score_above_reference_filtered"],
                supported["top15_above_reference_filtered"],
            ],
            ["both_above_reference", "score_only_above", "top15_only_above"],
            default="neither_above_reference",
        )
    else:
        reference_score = float("nan")
        reference_top15 = float("nan")
        supported["score_above_reference_filtered"] = False
        supported["top15_above_reference_filtered"] = False
        supported["reference_quadrant"] = "reference_unavailable"

    transfer_columns = [
        "step6_rank_filtered",
        "step7_rank_filtered",
        "step8_rank_filtered",
        "top15_rank_filtered",
        "score_rank_filtered",
        "score_rank_percentile",
        "top15_rank_percentile",
        "score_rank_minus_top15_rank",
        "step6_rank_minus_step8_rank",
        "top_score_list_flag",
        "top15_list_flag",
        "top_step8_list_flag",
        "top_list_membership",
        "is_reference_variant",
        "score_above_reference_filtered",
        "top15_above_reference_filtered",
        "reference_quadrant",
    ]
    boolean_columns = {
        "top_score_list_flag",
        "top15_list_flag",
        "top_step8_list_flag",
        "is_reference_variant",
        "score_above_reference_filtered",
        "top15_above_reference_filtered",
    }
    string_columns = {"top_list_membership", "reference_quadrant"}
    for column in transfer_columns:
        if column in boolean_columns:
            output[column] = False
        elif column in string_columns:
            output[column] = pd.Series(pd.NA, index=output.index, dtype="object")
        else:
            output[column] = np.nan
        output.loc[supported.index, column] = supported[column]

    total_only = output.loc[output["total_count_filter_pass"] & finite_metrics]
    summary_rows: list[dict[str, object]] = []
    correlation_pairs = [
        (
            "step6_high15_vs_step7_expected_score",
            "step6_high15_probability",
            "step7_expected_bin_score",
        ),
        (
            "step6_high15_vs_step8_high15_relative",
            "step6_high15_probability",
            "step8_high15_relative_enrichment",
        ),
        (
            "step7_expected_score_vs_step8_high15_relative",
            "step7_expected_bin_score",
            "step8_high15_relative_enrichment",
        ),
    ]

    def add_correlation_rows(scope: str, table: pd.DataFrame) -> None:
        for comparison, left_column, right_column in correlation_pairs:
            for method in ["spearman", "pearson", "kendall"]:
                statistic, pvalue, n = safe_correlation(
                    table[left_column], table[right_column], method
                )
                summary_rows.append(
                    {
                        "section": "correlation",
                        "scope": scope,
                        "comparison": comparison,
                        "method": method,
                        "metric": f"{method}_{comparison}",
                        "value": statistic,
                        "pvalue": pvalue,
                        "n": n,
                    }
                )
        # Preserve the pre-v0.2.2 metric name.
        statistic, pvalue, n = safe_correlation(
            table["step7_expected_bin_score"],
            table["step6_high15_probability"],
            "spearman",
        )
        summary_rows.append(
            {
                "section": "correlation",
                "scope": scope,
                "comparison": "legacy_expected_score_vs_high15_probability",
                "method": "spearman",
                "metric": "spearman_expected_score_vs_high15_probability",
                "value": statistic,
                "pvalue": pvalue,
                "n": n,
            }
        )
        if "top15_vs_unsorted_log2_enrichment" in table.columns:
            for method in ["spearman", "pearson", "kendall"]:
                statistic, pvalue, n = safe_correlation(
                    table["step7_expected_bin_score"],
                    table["top15_vs_unsorted_log2_enrichment"],
                    method,
                )
                summary_rows.append(
                    {
                        "section": "correlation",
                        "scope": scope,
                        "comparison": "expected_score_vs_unsorted_secondary",
                        "method": method,
                        "metric": (
                            f"{method}_expected_score_vs_"
                            "top15_unsorted_log2_secondary"
                        ),
                        "value": statistic,
                        "pvalue": pvalue,
                        "n": n,
                    }
                )

    add_correlation_rows("total_count_filter_only", total_only)
    add_correlation_rows("robust_read_support", supported)

    overlap_pairs = [
        ("step6_vs_step7", "step6_rank_filtered", "step7_rank_filtered"),
        ("step6_vs_step8", "step6_rank_filtered", "step8_rank_filtered"),
        ("step7_vs_step8", "step7_rank_filtered", "step8_rank_filtered"),
    ]
    overlap_rows: list[dict[str, object]] = []
    for requested_n in sorted({20, top_n, 100}):
        pair_results: dict[str, dict[str, object]] = {}
        for pair, left_rank, right_rank in overlap_pairs:
            overlap = rank_overlap(
                supported, requested_n, left_rank, right_rank, pair
            )
            overlap_rows.append(overlap)
            pair_results[pair] = overlap
        actual_k = min(requested_n, len(supported))
        if actual_k:
            three_sets = [
                set(supported.nsmallest(actual_k, rank)["variant_id"])
                for rank in [
                    "step6_rank_filtered",
                    "step7_rank_filtered",
                    "step8_rank_filtered",
                ]
            ]
            three_way_count = len(set.intersection(*three_sets))
            three_way_percent = 100 * three_way_count / actual_k
        else:
            three_way_count = 0
            three_way_percent = float("nan")
        overlap_rows.append(
            {
                "pair": "step6_step7_step8_three_way",
                "requested_top_n": requested_n,
                "actual_top_n": actual_k,
                "overlap_count": three_way_count,
                "overlap_percent_of_each_list": three_way_percent,
                "jaccard": float("nan"),
            }
        )
        # Preserve legacy summary names for the informative Step6-vs-Step7 pair.
        legacy = pair_results["step6_vs_step7"]
        for metric in [
            "actual_top_n",
            "overlap_count",
            "overlap_percent_of_each_list",
            "jaccard",
        ]:
            summary_rows.append(
                {
                    "section": "rank_overlap",
                    "scope": "robust_read_support",
                    "comparison": "step6_vs_step7",
                    "method": "top_n_overlap",
                    "metric": f"top{requested_n}_{metric}",
                    "value": legacy[metric],
                    "pvalue": np.nan,
                    "n": len(supported),
                }
            )

    step8_formula_error = (
        supported["step8_high15_relative_enrichment"]
        - supported["step6_high15_probability"] / high15_fraction
    ).abs()
    step6_step8_rank_mismatches = int(
        (supported["step6_rank_filtered"] != supported["step8_rank_filtered"]).sum()
    )
    counts = {
        "variants_total": int(len(output)),
        "variants_total_count_filter_pass": int(output["total_count_filter_pass"].sum()),
        "variants_comparison_read_support_pass": int(len(supported)),
        "configured_min_total_count": int(min_total_count),
        "configured_min_unsorted_count": int(min_unsorted_count),
        "configured_min_high_bin_count": int(min_high_bin_count),
        "configured_top_n": int(top_n),
        "actual_top_n": int(actual_top_n),
        "high15_nominal_population_fraction": high15_fraction,
        "top_list_consensus_count": int(
            supported["top_list_membership"].eq("consensus").sum()
        ),
        "top_list_score_only_count": int(
            supported["top_list_membership"].eq("score_only").sum()
        ),
        "top_list_top15_only_count": int(
            supported["top_list_membership"].eq("top15_only").sum()
        ),
        "reference_score": reference_score,
        "reference_high15_probability": reference_top15,
        "maximum_abs_step8_formula_error": (
            float(step8_formula_error.max()) if len(step8_formula_error) else float("nan")
        ),
        "step6_step8_rank_mismatch_count": step6_step8_rank_mismatches,
    }
    for key, value in counts.items():
        summary_rows.append(
            {
                "section": "count_or_setting",
                "scope": "all",
                "comparison": "not_applicable",
                "method": "not_applicable",
                "metric": key,
                "value": value,
                "pvalue": np.nan,
                "n": len(supported),
            }
        )

    if len(reference) == 1:
        quadrant_counts = supported["reference_quadrant"].value_counts()
        for quadrant in [
            "both_above_reference",
            "score_only_above",
            "top15_only_above",
            "neither_above_reference",
        ]:
            summary_rows.append(
                {
                    "section": "reference_quadrant",
                    "scope": "robust_read_support",
                    "comparison": "step6_vs_step7",
                    "method": "count",
                    "metric": quadrant,
                    "value": int(quadrant_counts.get(quadrant, 0)),
                    "pvalue": np.nan,
                    "n": len(supported),
                }
            )

    summary_table = pd.DataFrame(summary_rows)
    overlap_table = pd.DataFrame(overlap_rows)

    def summary_correlation(comparison: str, method: str = "spearman") -> float:
        selected = summary_table[
            (summary_table["scope"] == "robust_read_support")
            & (summary_table["comparison"] == comparison)
            & (summary_table["method"] == method)
        ]
        return float(selected["value"].iloc[0]) if len(selected) else float("nan")

    primary_rho = summary_correlation("step6_high15_vs_step7_expected_score")
    step6_step8_rho = summary_correlation(
        "step6_high15_vs_step8_high15_relative"
    )
    selected_overlap = rank_overlap(
        supported,
        top_n,
        "step6_rank_filtered",
        "step7_rank_filtered",
        "step6_vs_step7",
    )
    overlap_percent = float(selected_overlap["overlap_percent_of_each_list"])
    if np.isfinite(primary_rho) and primary_rho >= 0.70 and overlap_percent >= 60:
        agreement = "strong"
    elif np.isfinite(primary_rho) and primary_rho >= 0.40:
        agreement = "moderate"
    else:
        agreement = "weak_or_endpoint_divergent"
    summary_json = {
        **counts,
        "population_fractions": fractions.tolist(),
        "primary_spearman_rho": primary_rho,
        "step6_vs_step7_spearman_rho": primary_rho,
        "step6_vs_step8_spearman_rho": step6_step8_rho,
        "top_n_overlap_percent": overlap_percent,
        "step6_vs_step7_top_n_overlap_percent": overlap_percent,
        "agreement_class": agreement,
        "step8_identity_check_pass": bool(
            step6_step8_rank_mismatches == 0
            and (len(step8_formula_error) == 0 or step8_formula_error.max() < 1e-12)
        ),
        "recommended_interpretation": (
            "Rank clones by Step 6 conditional High15 probability. Use Step 7 "
            "expected score as whole-profile support. Use Step 8 p/w to visualize "
            "which bins are enriched; its High15 scalar is a rescaling of Step 6, "
            "not an independent endpoint."
        ),
    }

    supported = supported.sort_values(
        ["top_list_membership", "step6_rank_filtered", "step7_rank_filtered"]
    )
    output = output.sort_values(
        ["comparison_read_support_pass", "step6_rank_filtered"],
        ascending=[False, True],
        na_position="last",
    )
    return output, supported, summary_table, overlap_table, summary_json


def main() -> int:
    args = parse_args()
    result = read_table(args.input)
    fractions = read_population_fractions(args.sample_map)
    output, supported, summary, overlap, summary_json = analyze_metric_comparison(
        result,
        min_total_count=args.min_total_count,
        min_unsorted_count=args.min_unsorted_count,
        min_high_bin_count=args.min_high_bin_count,
        top_n=args.top_n,
        population_fractions=fractions,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    output.to_csv(
        args.outdir / "metric_comparison_all.csv", index=False, encoding="utf-8-sig"
    )
    supported.to_csv(
        args.outdir / "metric_comparison_supported.csv",
        index=False,
        encoding="utf-8-sig",
    )
    summary.to_csv(
        args.outdir / "metric_comparison_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    overlap.to_csv(
        args.outdir / "metric_comparison_topn_overlap.csv",
        index=False,
        encoding="utf-8-sig",
    )

    # Explicit Step 6/7/8 files make the new outputs discoverable while the
    # legacy filenames above remain available to older notebooks and reports.
    supported.to_csv(
        args.outdir / "step6_step7_step8_metrics.csv",
        index=False,
        encoding="utf-8-sig",
    )
    summary.loc[summary["section"].eq("correlation")].to_csv(
        args.outdir / "step6_step7_step8_correlations.csv",
        index=False,
        encoding="utf-8-sig",
    )
    overlap.to_csv(
        args.outdir / "step6_step7_step8_topn_overlap.csv",
        index=False,
        encoding="utf-8-sig",
    )
    supported[supported["top_list_membership"] == "consensus"].to_csv(
        args.outdir / "top_candidates_consensus.csv",
        index=False,
        encoding="utf-8-sig",
    )
    supported[supported["top_list_membership"] == "score_only"].to_csv(
        args.outdir / "top_candidates_score_only.csv",
        index=False,
        encoding="utf-8-sig",
    )
    supported[supported["top_list_membership"] == "top15_only"].to_csv(
        args.outdir / "top_candidates_top15_only.csv",
        index=False,
        encoding="utf-8-sig",
    )
    json_text = json.dumps(summary_json, indent=2, ensure_ascii=False)
    (args.outdir / "metric_comparison_summary.json").write_text(
        json_text, encoding="utf-8"
    )
    (args.outdir / "step6_step7_step8_summary.json").write_text(
        json_text, encoding="utf-8"
    )
    print(f"Step 6/7/8 comparison written to: {args.outdir}")
    print(
        "Shared comparison universe: "
        f"{summary_json['variants_comparison_read_support_pass']:,} UTRs"
    )
    print(
        "Step 6 vs Step 7 Spearman rho: "
        f"{summary_json['step6_vs_step7_spearman_rho']:.4f}"
    )
    print(
        f"Step 6 vs Step 7 Top {summary_json['actual_top_n']} overlap: "
        f"{summary_json['top_list_consensus_count']:,} "
        f"({summary_json['step6_vs_step7_top_n_overlap_percent']:.2f}%)"
    )
    print(
        "Step 6 vs Step 8 identity: "
        f"rho={summary_json['step6_vs_step8_spearman_rho']:.4f}, "
        f"rank mismatches={summary_json['step6_step8_rank_mismatch_count']}"
    )
    print(f"Agreement class (Step 6 vs 7): {summary_json['agreement_class']}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
