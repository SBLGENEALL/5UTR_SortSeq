#!/usr/bin/env python3
"""Compare whole-distribution score with bin1+bin2 enrichment rankings."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, pearsonr, spearmanr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare expected-bin score and top-15% enrichment after a shared "
            "six-bin total-read filter."
        )
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--min-total-count", type=int, default=201)
    parser.add_argument("--min-unsorted-count", type=int, default=50)
    parser.add_argument("--min-high-bin-count", type=int, default=20)
    parser.add_argument("--top-n", type=int, default=50)
    return parser.parse_args()


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


def rank_overlap(table: pd.DataFrame, top_n: int) -> dict[str, object]:
    actual_n = min(top_n, len(table))
    if actual_n == 0:
        return {
            "requested_top_n": top_n,
            "actual_top_n": 0,
            "overlap_count": 0,
            "overlap_percent_of_each_list": float("nan"),
            "jaccard": float("nan"),
        }
    score_ids = set(table.nsmallest(actual_n, "score_rank_filtered")["variant_id"])
    top15_ids = set(table.nsmallest(actual_n, "top15_rank_filtered")["variant_id"])
    overlap = score_ids & top15_ids
    union = score_ids | top15_ids
    return {
        "requested_top_n": top_n,
        "actual_top_n": actual_n,
        "overlap_count": len(overlap),
        "overlap_percent_of_each_list": 100 * len(overlap) / actual_n,
        "jaccard": len(overlap) / len(union) if union else float("nan"),
    }


def analyze_metric_comparison(
    result: pd.DataFrame,
    min_total_count: int = 201,
    min_unsorted_count: int = 50,
    min_high_bin_count: int = 20,
    top_n: int = 50,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    dict[str, object],
]:
    probability_columns = [f"bin{x}_probability" for x in range(1, 7)]
    required = {
        "variant_id",
        "total_6bin_count",
        "unsorted_count",
        "high_bin_raw_count",
        "expected_bin_score",
        "high15_probability",
        "top15_vs_unsorted_log2_enrichment",
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
        "top15_vs_unsorted_log2_enrichment",
        *probability_columns,
    ]
    for column in numeric_columns:
        output[column] = pd.to_numeric(output[column], errors="coerce")

    finite_metrics = (
        np.isfinite(output["expected_bin_score"])
        & np.isfinite(output["top15_vs_unsorted_log2_enrichment"])
    )
    output["total_count_filter_pass"] = (
        output["total_6bin_count"] >= min_total_count
    )
    output["comparison_read_support_pass"] = (
        output["total_count_filter_pass"]
        & (output["unsorted_count"] >= min_unsorted_count)
        & (output["high_bin_raw_count"] >= min_high_bin_count)
        & finite_metrics
    )

    supported_mask = output["comparison_read_support_pass"]
    supported = output.loc[supported_mask].copy()
    supported["score_rank_filtered"] = supported["expected_bin_score"].rank(
        method="min", ascending=False
    )
    supported["top15_rank_filtered"] = supported[
        "top15_vs_unsorted_log2_enrichment"
    ].rank(method="min", ascending=False)
    supported["score_rank_percentile"] = supported["expected_bin_score"].rank(
        method="average", pct=True
    )
    supported["top15_rank_percentile"] = supported[
        "top15_vs_unsorted_log2_enrichment"
    ].rank(method="average", pct=True)
    supported["score_rank_minus_top15_rank"] = (
        supported["score_rank_filtered"] - supported["top15_rank_filtered"]
    )

    actual_top_n = min(top_n, len(supported))
    supported["top_score_list_flag"] = False
    supported["top15_list_flag"] = False
    if actual_top_n:
        score_top_index = supported.nsmallest(
            actual_top_n, "score_rank_filtered"
        ).index
        top15_top_index = supported.nsmallest(
            actual_top_n, "top15_rank_filtered"
        ).index
        supported.loc[score_top_index, "top_score_list_flag"] = True
        supported.loc[top15_top_index, "top15_list_flag"] = True
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
        reference_score = float(reference["expected_bin_score"].iloc[0])
        reference_top15 = float(
            reference["top15_vs_unsorted_log2_enrichment"].iloc[0]
        )
        supported["score_above_reference_filtered"] = (
            supported["expected_bin_score"] > reference_score
        )
        supported["top15_above_reference_filtered"] = (
            supported["top15_vs_unsorted_log2_enrichment"] > reference_top15
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

    boolean_output_columns = {
        "top_score_list_flag",
        "top15_list_flag",
        "is_reference_variant",
        "score_above_reference_filtered",
        "top15_above_reference_filtered",
    }
    string_output_columns = {"top_list_membership", "reference_quadrant"}
    for column in [
        "score_rank_filtered",
        "top15_rank_filtered",
        "score_rank_percentile",
        "top15_rank_percentile",
        "score_rank_minus_top15_rank",
        "top_score_list_flag",
        "top15_list_flag",
        "top_list_membership",
        "is_reference_variant",
        "score_above_reference_filtered",
        "top15_above_reference_filtered",
        "reference_quadrant",
    ]:
        if column in boolean_output_columns:
            output[column] = False
        elif column in string_output_columns:
            output[column] = pd.Series(pd.NA, index=output.index, dtype="object")
        else:
            output[column] = np.nan
        output.loc[supported.index, column] = supported[column]

    total_only = output.loc[output["total_count_filter_pass"] & finite_metrics]
    summary_rows: list[dict[str, object]] = []

    def add_correlation_rows(scope: str, table: pd.DataFrame) -> None:
        pairs = [
            (
                "expected_score_vs_top15_unsorted_log2",
                table["expected_bin_score"],
                table["top15_vs_unsorted_log2_enrichment"],
            ),
            (
                "expected_score_vs_high15_probability",
                table["expected_bin_score"],
                table["high15_probability"],
            ),
        ]
        for comparison, left, right in pairs:
            for method in ["spearman", "pearson", "kendall"]:
                statistic, pvalue, n = safe_correlation(left, right, method)
                summary_rows.append(
                    {
                        "section": "correlation",
                        "scope": scope,
                        "metric": f"{method}_{comparison}",
                        "value": statistic,
                        "pvalue": pvalue,
                        "n": n,
                    }
                )

    add_correlation_rows("total_count_filter_only", total_only)
    add_correlation_rows("robust_read_support", supported)

    overlap_rows: list[dict[str, object]] = []
    for requested_n in sorted({20, top_n, 100}):
        overlap = rank_overlap(supported, requested_n)
        overlap_rows.append(overlap)
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
                    "metric": f"top{requested_n}_{metric}",
                    "value": overlap[metric],
                    "pvalue": np.nan,
                    "n": len(supported),
                }
            )

    counts = {
        "variants_total": int(len(output)),
        "variants_total_count_filter_pass": int(
            output["total_count_filter_pass"].sum()
        ),
        "variants_comparison_read_support_pass": int(len(supported)),
        "configured_min_total_count": int(min_total_count),
        "configured_min_unsorted_count": int(min_unsorted_count),
        "configured_min_high_bin_count": int(min_high_bin_count),
        "configured_top_n": int(top_n),
        "actual_top_n": int(actual_top_n),
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
        "reference_top15_log2_enrichment": reference_top15,
    }
    for key, value in counts.items():
        summary_rows.append(
            {
                "section": "count_or_setting",
                "scope": "all",
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
                    "metric": quadrant,
                    "value": int(quadrant_counts.get(quadrant, 0)),
                    "pvalue": np.nan,
                    "n": len(supported),
                }
            )

    summary_table = pd.DataFrame(summary_rows)
    overlap_table = pd.DataFrame(overlap_rows)
    primary_rho_row = summary_table[
        (summary_table["scope"] == "robust_read_support")
        & (
            summary_table["metric"]
            == "spearman_expected_score_vs_top15_unsorted_log2"
        )
    ]
    primary_rho = (
        float(primary_rho_row["value"].iloc[0])
        if len(primary_rho_row)
        else float("nan")
    )
    selected_overlap = rank_overlap(supported, top_n)
    overlap_percent = float(selected_overlap["overlap_percent_of_each_list"])
    if np.isfinite(primary_rho) and primary_rho >= 0.70 and overlap_percent >= 60:
        agreement = "strong"
    elif np.isfinite(primary_rho) and primary_rho >= 0.40:
        agreement = "moderate"
    else:
        agreement = "weak_or_endpoint_divergent"
    summary_json = {
        **counts,
        "primary_spearman_rho": primary_rho,
        "top_n_overlap_percent": overlap_percent,
        "agreement_class": agreement,
        "recommended_interpretation": (
            "Use top15 enrichment as the primary endpoint when the goal is the "
            "extreme-high fraction; use expected score for the average full-bin "
            "position. Prioritize consensus candidates above reference in both, "
            "then review discordant profiles and low-tail probability."
        ),
    }

    supported = supported.sort_values(
        ["top_list_membership", "top15_rank_filtered", "score_rank_filtered"]
    )
    output = output.sort_values(
        ["comparison_read_support_pass", "top15_rank_filtered"],
        ascending=[False, True],
        na_position="last",
    )
    return output, supported, summary_table, overlap_table, summary_json


def main() -> int:
    args = parse_args()
    separator = "\t" if args.input.suffix.lower() in {".tsv", ".txt"} else ","
    result = pd.read_csv(args.input, sep=separator)
    output, supported, summary, overlap, summary_json = analyze_metric_comparison(
        result,
        min_total_count=args.min_total_count,
        min_unsorted_count=args.min_unsorted_count,
        min_high_bin_count=args.min_high_bin_count,
        top_n=args.top_n,
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
    (args.outdir / "metric_comparison_summary.json").write_text(
        json.dumps(summary_json, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"Metric comparison written to: {args.outdir}")
    print(
        "Shared comparison universe: "
        f"{summary_json['variants_comparison_read_support_pass']:,} UTRs"
    )
    print(f"Spearman rho: {summary_json['primary_spearman_rho']:.4f}")
    print(
        f"Top {summary_json['actual_top_n']} overlap: "
        f"{summary_json['top_list_consensus_count']:,} "
        f"({summary_json['top_n_overlap_percent']:.2f}%)"
    )
    print(f"Agreement class: {summary_json['agreement_class']}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
