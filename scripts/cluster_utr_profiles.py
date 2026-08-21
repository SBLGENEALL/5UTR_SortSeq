#!/usr/bin/env python3
"""Export and cluster sum-to-one normalized enrichment profiles.

For UTR i and bin b, the displayed profile is

    q_ib = f_ib / sum_k(f_ik) = R_ib / sum_k(R_ik)

where f is sequencing-depth-normalized frequency and R=p/w is relative
enrichment. q sums to one across the six bins but is deliberately not called
a cell probability; it is an equal-bin normalized enrichment shape.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage


PROBABILITY_COLUMNS = [f"bin{x}_probability" for x in range(1, 7)]
RELATIVE_COLUMNS = [f"bin{x}_relative_enrichment" for x in range(1, 7)]
NORMALIZED_COLUMNS = [
    f"bin{x}_normalized_enrichment_share" for x in range(1, 7)
]
DEFAULT_POPULATION_FRACTIONS = np.asarray(
    [0.05, 0.10, 0.15, 0.20, 0.30, 0.20], dtype=float
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Export every UTR's q=f/sum(f) profile and group read-supported "
            "UTRs by Hellinger-transformed Ward clustering."
        )
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--sample-map", type=Path)
    parser.add_argument("--clusters", type=int, default=8)
    parser.add_argument("--top-n", type=int, default=200)
    parser.add_argument(
        "--min-total-count",
        type=int,
        default=200,
        help=(
            "Minimum sum of bin1-bin6 raw reads for inclusion in profile "
            "clustering (default: 200). All UTRs are still exported."
        ),
    )
    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    separator = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=separator)


def read_population_fractions(sample_map: Path | None) -> np.ndarray:
    if sample_map is None:
        return DEFAULT_POPULATION_FRACTIONS.copy()
    mapping = read_table(sample_map)
    fraction_column = (
        "population_fraction_used"
        if "population_fraction_used" in mapping.columns
        else "population_fraction"
    )
    required = {"sample_type", "bin_number", fraction_column}
    missing = sorted(required.difference(mapping.columns))
    if missing:
        raise ValueError(f"sample map is missing columns: {missing}")
    bins = mapping[
        mapping["sample_type"].astype(str).str.lower().eq("bin")
    ].copy()
    bins["bin_number"] = pd.to_numeric(
        bins["bin_number"], errors="raise"
    ).astype(int)
    bins = bins.sort_values("bin_number")
    if bins["bin_number"].tolist() != [1, 2, 3, 4, 5, 6]:
        raise ValueError("sample map must contain bin_number 1 through 6 exactly once")
    fractions = pd.to_numeric(
        bins[fraction_column], errors="raise"
    ).to_numpy(dtype=float)
    if fractions.sum() > 1.5:
        fractions = fractions / 100.0
    if (fractions <= 0).any() or not np.isclose(
        fractions.sum(), 1.0, atol=0.02
    ):
        raise ValueError("population fractions must be positive and sum to 1 or 100")
    return fractions / fractions.sum()


def as_bool(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "t", "1"})


def prepare_profiles(
    result: pd.DataFrame,
    population_fractions: np.ndarray,
) -> tuple[pd.DataFrame, str]:
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
    if len(population_fractions) != 6 or (population_fractions <= 0).any():
        raise ValueError("six positive population fractions are required")

    table = result.copy()
    table["variant_id"] = table["variant_id"].astype(str)
    if table["variant_id"].duplicated().any():
        raise ValueError("variant_id values must be unique")

    numeric_columns = [
        "expected_bin_score",
        "high15_probability",
        "high15_final_rank",
        "high15_primary_rank",
        "high_bin_raw_count",
        "total_6bin_count",
        "unsorted_count",
        *PROBABILITY_COLUMNS,
        *RELATIVE_COLUMNS,
        *NORMALIZED_COLUMNS,
    ]
    for column in numeric_columns:
        if column in table.columns:
            table[column] = pd.to_numeric(table[column], errors="coerce")

    probability = table[PROBABILITY_COLUMNS].to_numpy(dtype=float)
    relative = probability / population_fractions[np.newaxis, :]
    for index, column in enumerate(RELATIVE_COLUMNS):
        if column not in table.columns:
            table[column] = relative[:, index]

    if all(column in table.columns for column in NORMALIZED_COLUMNS):
        normalized = table[NORMALIZED_COLUMNS].to_numpy(dtype=float)
        source = "existing q columns from analyze_sortseq.py"
    else:
        normalized = relative.copy()
        source = "computed as (p/w)/sum(p/w), equivalent to f/sum(f)"

    row_sum = np.nansum(normalized, axis=1)
    usable = np.isfinite(normalized).all(axis=1) & (normalized >= 0).all(axis=1)
    usable &= np.isfinite(row_sum) & (row_sum > 0)
    normalized_output = np.full_like(normalized, np.nan, dtype=float)
    normalized_output[usable] = normalized[usable] / row_sum[usable, np.newaxis]
    for index, column in enumerate(NORMALIZED_COLUMNS):
        table[column] = normalized_output[:, index]
    table["normalized_enrichment_share_sum"] = np.nansum(
        normalized_output, axis=1
    )
    table.loc[~usable, "normalized_enrichment_share_sum"] = np.nan
    table["equal_bin_high_share"] = table[NORMALIZED_COLUMNS[:2]].sum(
        axis=1, min_count=2
    )

    if "is_reference_variant" in table.columns:
        table["is_reference_variant"] = as_bool(table["is_reference_variant"])
    else:
        table["is_reference_variant"] = table["variant_id"].str.lower().isin(
            {"original", "orginal"}
        )
    if "recommended_for_cloning" in table.columns:
        table["recommended_for_cloning"] = as_bool(
            table["recommended_for_cloning"]
        )
    else:
        table["recommended_for_cloning"] = False
    if "candidate_tier" not in table.columns:
        table["candidate_tier"] = "not_available"
    return table, source


def build_profile_outputs(
    result: pd.DataFrame,
    requested_clusters: int = 8,
    min_total_count: int = 200,
    population_fractions: np.ndarray | None = None,
    top_n: int = 200,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    if requested_clusters < 2:
        raise ValueError("clusters must be at least 2")
    if min_total_count < 0:
        raise ValueError("min_total_count must be nonnegative")
    if top_n <= 0:
        raise ValueError("top_n must be positive")
    fractions = (
        DEFAULT_POPULATION_FRACTIONS.copy()
        if population_fractions is None
        else np.asarray(population_fractions, dtype=float)
    )
    fractions = fractions / fractions.sum()
    table, normalized_source = prepare_profiles(result, fractions)

    finite_profile = np.isfinite(table[NORMALIZED_COLUMNS]).all(axis=1)
    nonnegative = (table[NORMALIZED_COLUMNS] >= 0).all(axis=1)
    positive_sum = table[NORMALIZED_COLUMNS].sum(axis=1) > 0
    total_support = (
        np.isfinite(table["total_6bin_count"])
        & (table["total_6bin_count"] >= min_total_count)
    )
    table["profile_eligible"] = (
        total_support & finite_profile & nonnegative & positive_sum
    )
    included = table.loc[table["profile_eligible"]].copy()
    if len(included) < 2:
        raise ValueError("at least two eligible UTR profiles are required")

    normalized = included[NORMALIZED_COLUMNS].to_numpy(dtype=float)
    hellinger = np.sqrt(normalized)
    hierarchy = linkage(hellinger, method="ward", optimal_ordering=True)
    actual_requested = min(requested_clusters, len(included))
    raw_cluster = fcluster(hierarchy, t=actual_requested, criterion="maxclust")
    included["raw_profile_cluster"] = raw_cluster

    cluster_order = (
        included.groupby("raw_profile_cluster", sort=False)
        .agg(
            median_equal_bin_high_share=("equal_bin_high_share", "median"),
            median_expected_bin_score=("expected_bin_score", "median"),
        )
        .sort_values(
            ["median_equal_bin_high_share", "median_expected_bin_score"],
            ascending=[False, False],
        )
        .index.tolist()
    )
    relabel = {raw: number for number, raw in enumerate(cluster_order, start=1)}
    included["profile_cluster"] = included["raw_profile_cluster"].map(
        relabel
    ).astype(int)

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

    counts = included["profile_cluster"].value_counts().sort_index()
    included["profile_label"] = included["profile_cluster"].map(
        lambda number: f"Profile {number} (n={int(counts.loc[number])})"
    )

    summary_rows: list[dict[str, object]] = []
    for cluster_number, group in included.groupby("profile_cluster", sort=True):
        mean_normalized = group[NORMALIZED_COLUMNS].mean()
        median_normalized = group[NORMALIZED_COLUMNS].median()
        mean_probability = group[PROBABILITY_COLUMNS].mean()
        median_probability = group[PROBABILITY_COLUMNS].median()
        row: dict[str, object] = {
            "profile_cluster": int(cluster_number),
            "profile_label": str(group["profile_label"].iloc[0]),
            "utr_count": int(len(group)),
            "utr_percent": 100 * len(group) / len(included),
            "dominant_mean_normalized_enrichment_bin": str(
                mean_normalized.idxmax()
            ).replace("_normalized_enrichment_share", ""),
            "mean_equal_bin_high_share": float(
                group["equal_bin_high_share"].mean()
            ),
            "median_equal_bin_high_share": float(
                group["equal_bin_high_share"].median()
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
        for number in range(1, 7):
            normalized_column = f"bin{number}_normalized_enrichment_share"
            probability_column = f"bin{number}_probability"
            row[f"bin{number}_mean_normalized_enrichment_share"] = float(
                mean_normalized[normalized_column]
            )
            row[f"bin{number}_median_normalized_enrichment_share"] = float(
                median_normalized[normalized_column]
            )
            row[f"bin{number}_mean_probability"] = float(
                mean_probability[probability_column]
            )
            row[f"bin{number}_median_probability"] = float(
                median_probability[probability_column]
            )
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)

    assignment_columns = [
        "variant_id", "profile_cluster", "profile_label", "heatmap_order",
        "dendrogram_order", "profile_eligible", "is_reference_variant",
        "candidate_tier", "recommended_for_cloning", "high15_final_rank",
        "high15_probability", "equal_bin_high_share", "expected_bin_score",
        "total_6bin_count", "unsorted_count", *PROBABILITY_COLUMNS,
        *RELATIVE_COLUMNS, *NORMALIZED_COLUMNS,
        "normalized_enrichment_share_sum",
    ]
    assignments = included[
        [column for column in assignment_columns if column in included.columns]
    ].sort_values("heatmap_order")

    cluster_metadata = included[
        ["variant_id", "profile_cluster", "profile_label", "heatmap_order"]
    ]
    all_profiles = table.drop(
        columns=[
            column
            for column in ["profile_cluster", "profile_label", "heatmap_order"]
            if column in table.columns
        ]
    ).merge(cluster_metadata, on="variant_id", how="left")
    all_profiles["all_utr_result_order"] = np.arange(1, len(all_profiles) + 1)
    all_columns = [
        "variant_id", "all_utr_result_order", "profile_eligible",
        "profile_cluster", "profile_label", "heatmap_order",
        "is_reference_variant", "high15_final_rank", "high15_primary_rank",
        "high15_probability", "equal_bin_high_share", "expected_bin_score",
        "candidate_tier", "recommended_for_cloning", "total_6bin_count",
        "high_bin_raw_count", "unsorted_count", *PROBABILITY_COLUMNS,
        *RELATIVE_COLUMNS, *NORMALIZED_COLUMNS,
        "normalized_enrichment_share_sum",
    ]
    all_profiles = all_profiles[
        [column for column in all_columns if column in all_profiles.columns]
    ]

    if "top15_read_support_pass" in table.columns:
        primary_support = as_bool(table["top15_read_support_pass"])
    else:
        primary_support = table["profile_eligible"]
    if "high15_final_rank" in table.columns and np.isfinite(
        table["high15_final_rank"]
    ).any():
        rank_column = "high15_final_rank"
        candidate_pool = table.loc[
            primary_support
            & finite_profile
            & np.isfinite(table[rank_column])
            & ~table["is_reference_variant"]
        ].sort_values([rank_column, "variant_id"])
    else:
        rank_column = "high15_probability"
        candidate_pool = table.loc[
            primary_support
            & finite_profile
            & np.isfinite(table[rank_column])
            & ~table["is_reference_variant"]
        ].sort_values([rank_column, "variant_id"], ascending=[False, True])
    top_candidates = candidate_pool.head(top_n).copy()
    top_candidates["top_candidate_selected"] = True
    top_candidates["reference_added_for_plot"] = False
    reference_rows = table.loc[table["is_reference_variant"]].copy()
    if len(reference_rows) > 0:
        reference_rows["top_candidate_selected"] = False
        reference_rows["reference_added_for_plot"] = True
        top_profiles = pd.concat([top_candidates, reference_rows], ignore_index=True)
    else:
        top_profiles = top_candidates.copy()
    top_profiles["top_profile_order"] = np.arange(1, len(top_profiles) + 1)
    top_columns = [
        "variant_id", "top_profile_order", "top_candidate_selected",
        "reference_added_for_plot", "is_reference_variant",
        "high15_final_rank", "high15_primary_rank", "high15_probability",
        "equal_bin_high_share", "expected_bin_score", "candidate_tier",
        "recommended_for_cloning", "total_6bin_count", "high_bin_raw_count",
        "unsorted_count", *PROBABILITY_COLUMNS, *RELATIVE_COLUMNS,
        *NORMALIZED_COLUMNS, "normalized_enrichment_share_sum",
    ]
    top_profiles = top_profiles[
        [column for column in top_columns if column in top_profiles.columns]
    ]

    reference_summary = reference_rows.iloc[0] if len(reference_rows) else None
    reference_cluster = None
    if len(reference_rows) and reference_rows.index[0] in included.index:
        reference_cluster = int(
            included.loc[reference_rows.index[0], "profile_cluster"]
        )
    manifest = {
        "variants_input": int(len(table)),
        "variants_exported_all": int(len(all_profiles)),
        "variants_included_in_clustering": int(len(included)),
        "variants_excluded_from_clustering": int(len(table) - len(included)),
        "variants_included": int(len(included)),
        "variants_excluded": int(len(table) - len(included)),
        "eligibility_column": "total_6bin_count",
        "minimum_total_6bin_count": int(min_total_count),
        "eligibility_rule": f"total_6bin_count >= {min_total_count}",
        "high_bin_support_filter_applied_to_clustering": False,
        "unsorted_filter_applied_to_clustering": False,
        "requested_clusters": int(requested_clusters),
        "actual_clusters": int(summary["profile_cluster"].nunique()),
        "profile_definition": "q_ib=f_ib/sum_k(f_ik)=R_ib/sum_k(R_ik)",
        "profile_values_sum_to_one": True,
        "profile_value_interpretation": (
            "normalized relative-enrichment shape, not P(bin|UTR)"
        ),
        "normalized_profile_source": normalized_source,
        "population_fractions": [float(value) for value in fractions],
        "distance": "Hellinger distance on six-bin q vectors",
        "linkage": "Ward",
        "cluster_label_order": (
            "median equal_bin_high_share descending, then expected score"
        ),
        "top_candidates_requested": int(top_n),
        "top_candidates_selected": int(len(top_candidates)),
        "top_candidate_rank_source": rank_column,
        "reference_variant_ids": reference_rows["variant_id"].tolist(),
        "reference_high15_final_rank": (
            None
            if reference_summary is None
            or "high15_final_rank" not in reference_summary
            or not np.isfinite(reference_summary["high15_final_rank"])
            else float(reference_summary["high15_final_rank"])
        ),
        "reference_profile_cluster": reference_cluster,
        "maximum_normalized_profile_sum_error": float(
            np.nanmax(
                np.abs(
                    table["normalized_enrichment_share_sum"].to_numpy(dtype=float)
                    - 1.0
                )
            )
        ),
    }
    return assignments, summary, all_profiles, top_profiles, manifest


def cluster_profiles(
    result: pd.DataFrame,
    requested_clusters: int = 8,
    min_total_count: int = 200,
    population_fractions: np.ndarray | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    """Backward-compatible three-output wrapper used by older callers."""
    assignments, summary, _, _, manifest = build_profile_outputs(
        result,
        requested_clusters=requested_clusters,
        min_total_count=min_total_count,
        population_fractions=population_fractions,
    )
    return assignments, summary, manifest


def main() -> int:
    args = parse_args()
    result = read_table(args.input)
    fractions = read_population_fractions(args.sample_map)
    assignments, summary, all_profiles, top_profiles, manifest = (
        build_profile_outputs(
            result,
            requested_clusters=args.clusters,
            min_total_count=args.min_total_count,
            population_fractions=fractions,
            top_n=args.top_n,
        )
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
    all_profiles.to_csv(
        args.outdir / "all_utr_normalized_enrichment_profiles.csv",
        index=False,
        encoding="utf-8-sig",
    )
    top_profiles.to_csv(
        args.outdir / "top_normalized_enrichment_profiles.csv",
        index=False,
        encoding="utf-8-sig",
    )
    top_profiles.to_csv(
        args.outdir / f"top{args.top_n}_normalized_enrichment_profiles.csv",
        index=False,
        encoding="utf-8-sig",
    )
    (args.outdir / "all_utr_profile_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"Normalized enrichment profiles written to: {args.outdir}")
    print(
        f"All UTRs exported: {manifest['variants_exported_all']:,}; "
        f"clustered: {manifest['variants_included_in_clustering']:,}; "
        f"clusters: {manifest['actual_clusters']}"
    )
    print(
        f"Top candidates for plots: {manifest['top_candidates_selected']:,} "
        f"plus reference when detected"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
