#!/usr/bin/env python3
"""Cluster read-supported UTRs by their complete six-bin probability profiles."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage


PROBABILITY_COLUMNS = [f"bin{x}_probability" for x in range(1, 7)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Group read-supported UTRs with similar six-bin conditional "
            "probability profiles using Hellinger-transformed Ward clustering."
        )
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--clusters", type=int, default=8)
    parser.add_argument(
        "--min-total-count",
        type=int,
        default=200,
        help=(
            "Minimum sum of bin1-bin6 raw reads for inclusion in the "
            "whole-library profile overview (default: 200)."
        ),
    )
    return parser.parse_args()


def as_bool(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "t", "1"})


def cluster_profiles(
    result: pd.DataFrame,
    requested_clusters: int = 8,
    min_total_count: int = 200,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    required = {
        "variant_id",
        "expected_bin_score",
        "high15_probability",
        "total_6bin_count",
        *PROBABILITY_COLUMNS,
    }
    missing = sorted(required.difference(result.columns))
    if missing:
        raise ValueError(f"input is missing columns: {missing}")
    if requested_clusters < 2:
        raise ValueError("clusters must be at least 2")
    if min_total_count < 0:
        raise ValueError("min_total_count must be nonnegative")

    table = result.copy()
    table["variant_id"] = table["variant_id"].astype(str)
    if table["variant_id"].duplicated().any():
        raise ValueError("variant_id values must be unique")
    numeric_columns = [
        "expected_bin_score",
        "high15_probability",
        "high15_final_rank",
        "total_6bin_count",
        "unsorted_count",
        *PROBABILITY_COLUMNS,
    ]
    for column in numeric_columns:
        if column in table.columns:
            table[column] = pd.to_numeric(table[column], errors="coerce")

    # The whole-library profile overview must contain low-, middle-, and
    # high-expression shapes. Therefore it uses only the total six-bin read
    # floor and deliberately does not inherit the High15-specific bin1+bin2
    # support filter or the unsorted QC flag.
    eligibility_column = "total_6bin_count"
    eligible = (
        np.isfinite(table["total_6bin_count"])
        & (table["total_6bin_count"] >= min_total_count)
    )
    finite = (
        np.isfinite(table[PROBABILITY_COLUMNS]).all(axis=1)
        & np.isfinite(table["expected_bin_score"])
        & np.isfinite(table["high15_probability"])
    )
    nonnegative = (table[PROBABILITY_COLUMNS] >= 0).all(axis=1)
    positive_sum = table[PROBABILITY_COLUMNS].sum(axis=1) > 0
    included = table.loc[eligible & finite & nonnegative & positive_sum].copy()
    if len(included) < 2:
        raise ValueError("at least two eligible UTR profiles are required")

    probability = included[PROBABILITY_COLUMNS].to_numpy(dtype=float)
    original_row_sums = probability.sum(axis=1)
    probability = probability / original_row_sums[:, np.newaxis]
    included.loc[:, PROBABILITY_COLUMNS] = probability

    # Square-rooting a compositional probability vector makes ordinary
    # Euclidean distance equal to Hellinger distance up to a constant factor.
    # Ward clustering then groups profile shape without using UTR abundance.
    hellinger = np.sqrt(probability)
    hierarchy = linkage(hellinger, method="ward", optimal_ordering=True)
    actual_requested = min(requested_clusters, len(included))
    raw_cluster = fcluster(hierarchy, t=actual_requested, criterion="maxclust")
    included["raw_profile_cluster"] = raw_cluster

    cluster_order = (
        included.groupby("raw_profile_cluster", sort=False)
        .agg(
            median_expected_bin_score=("expected_bin_score", "median"),
            median_high15_probability=("high15_probability", "median"),
        )
        .sort_values(
            ["median_expected_bin_score", "median_high15_probability"],
            ascending=[False, False],
        )
        .index.tolist()
    )
    relabel = {raw: number for number, raw in enumerate(cluster_order, start=1)}
    included["profile_cluster"] = included["raw_profile_cluster"].map(relabel).astype(int)

    leaf_position = np.empty(len(included), dtype=int)
    leaf_position[leaves_list(hierarchy)] = np.arange(1, len(included) + 1)
    included["dendrogram_order"] = leaf_position
    heatmap_index = included.sort_values(
        ["profile_cluster", "dendrogram_order"]
    ).index
    heatmap_order = pd.Series(
        np.arange(1, len(included) + 1), index=heatmap_index, dtype=int
    )
    included["heatmap_order"] = heatmap_order.reindex(included.index).to_numpy()

    if "is_reference_variant" in included.columns:
        included["is_reference_variant"] = as_bool(included["is_reference_variant"])
    else:
        included["is_reference_variant"] = included["variant_id"].str.lower().isin(
            {"original", "orginal"}
        )
    if "recommended_for_cloning" in included.columns:
        included["recommended_for_cloning"] = as_bool(
            included["recommended_for_cloning"]
        )
    else:
        included["recommended_for_cloning"] = False
    if "candidate_tier" not in included.columns:
        included["candidate_tier"] = "not_available"

    counts = included["profile_cluster"].value_counts().sort_index()
    included["profile_label"] = included["profile_cluster"].map(
        lambda number: f"Profile {number} (n={int(counts.loc[number])})"
    )

    summary_rows: list[dict[str, object]] = []
    for cluster_number, group in included.groupby("profile_cluster", sort=True):
        mean_probability = group[PROBABILITY_COLUMNS].mean()
        median_probability = group[PROBABILITY_COLUMNS].median()
        row: dict[str, object] = {
            "profile_cluster": int(cluster_number),
            "profile_label": str(group["profile_label"].iloc[0]),
            "utr_count": int(len(group)),
            "utr_percent": 100 * len(group) / len(included),
            "dominant_mean_bin": str(mean_probability.idxmax()).replace(
                "_probability", ""
            ),
            "mean_high15_probability": float(group["high15_probability"].mean()),
            "median_high15_probability": float(
                group["high15_probability"].median()
            ),
            "mean_expected_bin_score": float(group["expected_bin_score"].mean()),
            "median_expected_bin_score": float(
                group["expected_bin_score"].median()
            ),
            "reference_present": bool(group["is_reference_variant"].any()),
            "recommended_cloning_count": int(
                group["recommended_for_cloning"].sum()
            ),
        }
        for number, column in enumerate(PROBABILITY_COLUMNS, start=1):
            row[f"bin{number}_mean_probability"] = float(mean_probability[column])
            row[f"bin{number}_median_probability"] = float(
                median_probability[column]
            )
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)

    preferred_columns = [
        "variant_id",
        "profile_cluster",
        "profile_label",
        "heatmap_order",
        "dendrogram_order",
        "is_reference_variant",
        "candidate_tier",
        "recommended_for_cloning",
        "high15_final_rank",
        "high15_probability",
        "expected_bin_score",
        "total_6bin_count",
        "unsorted_count",
        *PROBABILITY_COLUMNS,
    ]
    assignments = included[
        [column for column in preferred_columns if column in included.columns]
    ].sort_values("heatmap_order")
    manifest = {
        "variants_input": int(len(table)),
        "variants_included": int(len(included)),
        "variants_excluded": int(len(table) - len(included)),
        "eligibility_column": eligibility_column,
        "minimum_total_6bin_count": int(min_total_count),
        "eligibility_rule": f"total_6bin_count >= {min_total_count}",
        "high_bin_support_filter_applied": False,
        "unsorted_filter_applied": False,
        "requested_clusters": int(requested_clusters),
        "actual_clusters": int(summary["profile_cluster"].nunique()),
        "distance": "Hellinger distance on six-bin probability vectors",
        "linkage": "Ward",
        "cluster_label_order": "median expected_bin_score descending",
        "maximum_original_probability_sum_error": float(
            np.max(np.abs(original_row_sums - 1.0))
        ),
    }
    return assignments, summary, manifest


def main() -> int:
    args = parse_args()
    separator = "\t" if args.input.suffix.lower() in {".tsv", ".txt"} else ","
    result = pd.read_csv(args.input, sep=separator)
    assignments, summary, manifest = cluster_profiles(
        result,
        requested_clusters=args.clusters,
        min_total_count=args.min_total_count,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    assignments.to_csv(
        args.outdir / "all_utr_profile_assignments.csv",
        index=False,
        encoding="utf-8-sig",
    )
    summary.to_csv(
        args.outdir / "all_utr_profile_cluster_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    (args.outdir / "all_utr_profile_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"All-UTR profile clustering written to: {args.outdir}")
    print(
        f"Included UTRs: {manifest['variants_included']:,}; "
        f"clusters: {manifest['actual_clusters']}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
