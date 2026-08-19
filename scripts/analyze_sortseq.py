#!/usr/bin/env python3
"""Easy analysis of one whole-population unsorted sample and six FACS bins.

Expected bin direction is fixed for clarity:
    bin1 = highest mCherry (5%)
    bin2 = next-highest (10%)
    ...
    bin6 = lowest mCherry-positive bin (20%)

For each UTR, the estimated fraction of its gated cells in bin b is:

    P(bin=b | UTR, gate) proportional to
        population_fraction[b] * within_sample_frequency[UTR, b]

This corrects both unequal NGS depth and unequal FACS-bin widths.
"""

from __future__ import annotations

import argparse
import html
import json
import math
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr


REQUIRED_SAMPLE_COLUMNS = {"sample_id", "sample_type", "count_file"}
ALLOWED_SAMPLE_TYPES = {"unsorted", "bin", "plasmid"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Analyze one whole-population unsorted sample plus six high-to-low "
            "mCherry FACS bins."
        )
    )
    parser.add_argument("--variants", required=True, type=Path)
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--id-column", default="variant_id")
    parser.add_argument("--pseudocount", type=float, default=0.5)
    parser.add_argument("--min-unsorted-count", type=int, default=50)
    parser.add_argument("--min-total-bin-count", type=int, default=100)
    parser.add_argument("--top-hit-min-unsorted-count", type=int, default=50)
    parser.add_argument("--top-hit-min-total-bin-count", type=int, default=200)
    parser.add_argument("--top-hit-min-high-bin-count", type=int, default=20)
    parser.add_argument("--top-hit-min-enrichment", type=float, default=1.0)
    parser.add_argument("--strict-min-unsorted-count", type=int, default=1000)
    parser.add_argument("--strict-min-total-bin-count", type=int, default=5000)
    parser.add_argument("--strict-relative-median-fraction", type=float, default=0.10)
    parser.add_argument("--strict-min-high-bin-count", type=int, default=200)
    parser.add_argument("--strict-min-detected-bins", type=int, default=3)
    parser.add_argument("--jackpot-max-bin-probability", type=float, default=0.85)
    parser.add_argument("--jackpot-max-detected-bins", type=int, default=2)
    parser.add_argument("--overall-gate-fraction", type=float, default=0.90)
    parser.add_argument("--top-n", type=int, default=50)
    parser.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=1000,
        help=(
            "Multinomial read-resampling replicates for technical High15 rank "
            "stability. This does not estimate biological or PCR uncertainty."
        ),
    )
    parser.add_argument("--bootstrap-seed", type=int, default=20260819)
    parser.add_argument("--bootstrap-lower-quantile", type=float, default=0.10)
    parser.add_argument("--bootstrap-top-n", type=int, default=50)
    parser.add_argument(
        "--bootstrap-min-probability-above-reference", type=float, default=0.90
    )
    parser.add_argument(
        "--reference-variant-id",
        default="auto",
        help=(
            "Reference/control variant ID. 'auto' recognizes original/orginal "
            "case-insensitively; 'none' disables reference comparison."
        ),
    )
    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    separator = "," if path.suffix.lower() == ".csv" else "\t"
    return pd.read_csv(path, sep=separator)


def gini(values: pd.Series) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0 or x.sum() <= 0:
        return float("nan")
    x = np.sort(np.maximum(x, 0))
    n = len(x)
    return float((2 * np.sum((np.arange(n) + 1) * x) / (n * x.sum())) - (n + 1) / n)


def load_variants(path: Path, id_column: str) -> pd.DataFrame:
    variants = read_table(path)
    if id_column not in variants.columns:
        raise ValueError(f"variants file must contain {id_column!r}")
    variants = variants.copy()
    variants[id_column] = variants[id_column].astype(str)
    if variants[id_column].duplicated().any():
        duplicated = variants.loc[variants[id_column].duplicated(), id_column].tolist()
        raise ValueError(f"variant IDs must be unique; duplicates include {duplicated[:5]}")
    return variants


def resolve_reference_variant_id(
    variants: pd.DataFrame,
    id_column: str,
    requested: str | None,
) -> str | None:
    """Resolve an explicit ID or the common original/orginal control aliases."""
    value = "auto" if requested is None else str(requested).strip()
    if value.lower() in {"", "none", "off", "false"}:
        return None
    identifiers = variants[id_column].astype(str)
    normalized = identifiers.str.strip().str.lower()
    aliases = {"original", "orginal"} if value.lower() == "auto" else {value.lower()}
    direct_mask = normalized.isin(aliases)
    original_ids_mask = pd.Series(False, index=variants.index)
    if "original_variant_ids" in variants.columns:
        original_ids_mask = variants["original_variant_ids"].fillna("").astype(str).map(
            lambda item: any(part.strip().lower() in aliases for part in item.split("|"))
        )
    matches = variants.loc[direct_mask | original_ids_mask, id_column]
    values = matches.astype(str).tolist()
    if len(values) == 0:
        if value.lower() == "auto":
            return None
        raise ValueError(f"reference variant ID not found: {value!r}")
    if len(values) > 1:
        raise ValueError(f"reference variant ID is ambiguous: {values}")
    return values[0]


def reference_display_label(
    variants: pd.DataFrame,
    id_column: str,
    reference_variant_id: str | None,
) -> str | None:
    if reference_variant_id is None:
        return None
    row = variants.loc[variants[id_column].astype(str).eq(reference_variant_id)].iloc[0]
    if "original_variant_ids" in variants.columns:
        identifiers = str(row["original_variant_ids"]).split("|")
        original_alias = next(
            (item.strip() for item in identifiers if item.strip().lower() in {"original", "orginal"}),
            None,
        )
        if original_alias:
            return original_alias
    return reference_variant_id


def load_samples(path: Path) -> pd.DataFrame:
    samples = read_table(path)
    missing = REQUIRED_SAMPLE_COLUMNS.difference(samples.columns)
    if missing:
        raise ValueError(f"samples file is missing columns: {sorted(missing)}")
    samples = samples.copy()
    samples["sample_id"] = samples["sample_id"].astype(str)
    samples["sample_type"] = samples["sample_type"].astype(str).str.lower()
    invalid = sorted(set(samples["sample_type"]) - ALLOWED_SAMPLE_TYPES)
    if invalid:
        raise ValueError(f"invalid sample_type values: {invalid}")
    if samples["sample_id"].duplicated().any():
        raise ValueError("sample_id values must be unique")

    for column in [
        "bin_number",
        "population_fraction",
        "cells_collected",
        "representative_log10_mfi",
    ]:
        if column not in samples.columns:
            samples[column] = np.nan
        samples[column] = pd.to_numeric(samples[column], errors="coerce")
    if "id_column" not in samples.columns:
        samples["id_column"] = "variant_id"
    if "count_column" not in samples.columns:
        samples["count_column"] = "count"
    samples["id_column"] = samples["id_column"].fillna("variant_id").astype(str)
    samples["count_column"] = samples["count_column"].fillna("count").astype(str)

    unsorted = samples[samples["sample_type"].eq("unsorted")]
    bins = samples[samples["sample_type"].eq("bin")]
    plasmid = samples[samples["sample_type"].eq("plasmid")]
    if len(unsorted) != 1:
        raise ValueError(f"exactly one unsorted sample is required; found {len(unsorted)}")
    if len(bins) != 6:
        raise ValueError(f"exactly six bin samples are required; found {len(bins)}")
    if len(plasmid) > 1:
        raise ValueError("at most one plasmid sample is allowed")
    if bins["bin_number"].isna().any():
        raise ValueError("bin_number is required for all six bin rows")
    bins_numbers = sorted(bins["bin_number"].astype(int).tolist())
    if bins_numbers != [1, 2, 3, 4, 5, 6]:
        raise ValueError("bin_number must contain 1,2,3,4,5,6 exactly once")
    if bins["population_fraction"].isna().any():
        raise ValueError("population_fraction is required for all six bins")

    fraction_index = bins.index
    # The gate-design fractions describe the proportion of the target-gated
    # population represented by each simultaneous sort bin.  Collected-event
    # counts can differ because of recovery/yield and are therefore retained as
    # QC metadata, not used to silently redefine the phenotype bins.
    fractions = bins["population_fraction"].astype(float).to_numpy()
    if fractions.sum() > 1.5:
        fractions = fractions / 100.0
    fraction_source = "population_fraction"

    samples["cells_collected_fraction"] = np.nan
    collected = bins["cells_collected"].astype(float)
    if collected.notna().all() and (collected > 0).all():
        samples.loc[fraction_index, "cells_collected_fraction"] = (
            collected.to_numpy() / collected.sum()
        )
    if (fractions <= 0).any() or not np.isclose(fractions.sum(), 1.0, atol=0.02):
        raise ValueError(
            f"population_fraction values must be positive and sum to 1 (or 100); "
            f"observed sum={fractions.sum():.4f}"
        )
    samples.loc[fraction_index, "population_fraction"] = fractions / fractions.sum()
    samples["population_fraction_source"] = "not_applicable"
    samples.loc[fraction_index, "population_fraction_source"] = fraction_source
    return samples


def technical_high15_bootstrap(
    raw_bins: pd.DataFrame,
    fractions: pd.Series,
    reference_variant_id: str | None,
    read_support: pd.Series,
    replicates: int,
    seed: int,
    lower_quantile: float,
    top_n: int,
) -> pd.DataFrame:
    """Estimate read-sampling stability of the conditional High15 endpoint.

    Each bin's complete observed variant composition is resampled with a
    multinomial distribution at the bin's original read depth.  This propagates
    finite NGS read sampling through depth normalization, bin-size correction,
    and within-UTR normalization.  It intentionally does *not* model biological
    replication, sorted-cell sampling, PCR jackpotting, or other wet-lab bias.
    """
    if replicates < 0:
        raise ValueError("bootstrap_replicates cannot be negative")
    if not 0 < lower_quantile < 0.5:
        raise ValueError("bootstrap_lower_quantile must be between 0 and 0.5")
    if top_n <= 0:
        raise ValueError("bootstrap_top_n must be positive")

    output = pd.DataFrame(index=raw_bins.index)
    statistic_columns = [
        "high15_bootstrap_median",
        "high15_bootstrap_lower",
        "high15_bootstrap_upper",
        "high15_log2_fold_bootstrap_median",
        "high15_log2_fold_bootstrap_lower",
        "high15_log2_fold_bootstrap_upper",
        "high15_bootstrap_probability_above_comparator",
        "high15_bootstrap_top_n_frequency",
    ]
    if replicates == 0:
        for column in statistic_columns:
            output[column] = np.nan
        return output

    depths = raw_bins.sum(axis=0).astype(np.int64)
    if (depths <= 0).any():
        raise ValueError("cannot bootstrap an empty bin")
    fraction_values = fractions.astype(float).to_numpy()
    rng = np.random.default_rng(seed)
    n_variants = len(raw_bins)
    high_mass = np.zeros((replicates, n_variants), dtype=np.float64)
    total_mass = np.zeros_like(high_mass)

    for bin_position, sample_id in enumerate(raw_bins.columns):
        depth = int(depths[sample_id])
        probabilities = raw_bins[sample_id].astype(float).to_numpy() / depth
        # Guard against harmless floating-point drift before multinomial input.
        probabilities = probabilities / probabilities.sum()
        resampled = rng.multinomial(depth, probabilities, size=replicates)
        mass = (resampled / depth) * fraction_values[bin_position]
        total_mass += mass
        if bin_position < 2:
            high_mass += mass

    with np.errstate(divide="ignore", invalid="ignore"):
        high15 = np.divide(
            high_mass,
            total_mass,
            out=np.full_like(high_mass, np.nan),
            where=total_mass > 0,
        )

    if reference_variant_id is not None:
        if reference_variant_id not in raw_bins.index:
            raise ValueError(
                f"bootstrap reference variant not found: {reference_variant_id!r}"
            )
        reference_position = raw_bins.index.get_loc(reference_variant_id)
        comparator = high15[:, reference_position][:, np.newaxis]
    else:
        comparator = np.full((replicates, 1), float(fraction_values[:2].sum()))

    epsilon = 1e-12
    with np.errstate(divide="ignore", invalid="ignore"):
        log2_fold = np.log2((high15 + epsilon) / (comparator + epsilon))

    upper_quantile = 1.0 - lower_quantile
    output["high15_bootstrap_median"] = np.nanmedian(high15, axis=0)
    output["high15_bootstrap_lower"] = np.nanquantile(
        high15, lower_quantile, axis=0
    )
    output["high15_bootstrap_upper"] = np.nanquantile(
        high15, upper_quantile, axis=0
    )
    output["high15_log2_fold_bootstrap_median"] = np.nanmedian(log2_fold, axis=0)
    output["high15_log2_fold_bootstrap_lower"] = np.nanquantile(
        log2_fold, lower_quantile, axis=0
    )
    output["high15_log2_fold_bootstrap_upper"] = np.nanquantile(
        log2_fold, upper_quantile, axis=0
    )
    output["high15_bootstrap_probability_above_comparator"] = np.nanmean(
        log2_fold > 0, axis=0
    )

    support_positions = np.flatnonzero(read_support.reindex(raw_bins.index).to_numpy())
    inclusion_count = np.zeros(n_variants, dtype=np.int64)
    actual_top_n = min(top_n, len(support_positions))
    if actual_top_n > 0:
        supported_values = high15[:, support_positions]
        if actual_top_n == len(support_positions):
            top_positions = np.broadcast_to(
                support_positions, (replicates, len(support_positions))
            )
        else:
            local_top = np.argpartition(
                -supported_values, actual_top_n - 1, axis=1
            )[:, :actual_top_n]
            top_positions = support_positions[local_top]
        np.add.at(inclusion_count, top_positions.ravel(), 1)
    output["high15_bootstrap_top_n_frequency"] = inclusion_count / replicates
    return output


def load_count_matrix(
    samples: pd.DataFrame,
    variants: pd.DataFrame,
    id_column: str,
    samples_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    variant_index = pd.Index(variants[id_column].astype(str), name=id_column)
    matrix = pd.DataFrame(index=variant_index)
    qc_rows: list[dict[str, object]] = []

    for _, sample in samples.iterrows():
        sample_id = str(sample["sample_id"])
        count_path = Path(str(sample["count_file"]))
        if not count_path.is_absolute():
            count_path = (samples_path.parent / count_path).resolve()
        if not count_path.exists():
            raise FileNotFoundError(f"{sample_id}: count file not found: {count_path}")
        table = read_table(count_path)
        sample_id_column = str(sample["id_column"])
        count_column = str(sample["count_column"])
        required = {sample_id_column, count_column}
        if not required.issubset(table.columns):
            raise ValueError(f"{sample_id}: count file must contain {sorted(required)}")
        table = table[[sample_id_column, count_column]].copy()
        table[sample_id_column] = table[sample_id_column].astype(str)
        table[count_column] = pd.to_numeric(table[count_column], errors="raise")
        if (table[count_column] < 0).any():
            raise ValueError(f"{sample_id}: counts cannot be negative")
        grouped = table.groupby(sample_id_column, sort=False)[count_column].sum()
        unknown_ids = grouped.index.difference(variant_index)
        known = grouped.reindex(variant_index, fill_value=0).astype(float)
        matrix[sample_id] = known.to_numpy()
        qc_rows.append(
            {
                "sample_id": sample_id,
                "sample_type": sample["sample_type"],
                "bin_number": sample["bin_number"],
                "population_fraction_used": sample["population_fraction"],
                "population_fraction_source": sample["population_fraction_source"],
                "cells_collected_reported": sample["cells_collected"],
                "cells_collected_fraction_qc": sample["cells_collected_fraction"],
                "total_reads_in_count_file": float(grouped.sum()),
                "assigned_reference_reads": float(known.sum()),
                "detected_reference_variants": int((known > 0).sum()),
                "variants_with_at_least_50_reads": int((known >= 50).sum()),
                "unknown_variant_ids": int(len(unknown_ids)),
                "unknown_reads": float(grouped.reindex(unknown_ids).sum())
                if len(unknown_ids)
                else 0.0,
                "count_gini": gini(known),
            }
        )
    return matrix, pd.DataFrame(qc_rows)


def percentile_tier(percentile: float) -> str:
    if not np.isfinite(percentile):
        return "low_coverage"
    if percentile >= 99:
        return "top_1pct"
    if percentile >= 95:
        return "top_5pct"
    if percentile >= 90:
        return "top_10pct"
    if percentile >= 75:
        return "upper_25pct"
    if percentile >= 25:
        return "middle_50pct"
    if percentile >= 10:
        return "lower_25pct"
    return "bottom_10pct"


def analyze(
    variants: pd.DataFrame,
    samples: pd.DataFrame,
    counts: pd.DataFrame,
    id_column: str,
    pseudocount: float,
    min_unsorted_count: int,
    min_total_bin_count: int,
    overall_gate_fraction: float,
    reference_variant_id: str | None = None,
    reference_label: str | None = None,
    strict_min_unsorted_count: int = 1000,
    strict_min_total_bin_count: int = 5000,
    strict_relative_median_fraction: float = 0.10,
    strict_min_high_bin_count: int = 200,
    strict_min_detected_bins: int = 3,
    jackpot_max_bin_probability: float = 0.85,
    jackpot_max_detected_bins: int = 2,
    top_hit_min_unsorted_count: int = 50,
    top_hit_min_total_bin_count: int = 200,
    top_hit_min_high_bin_count: int = 20,
    top_hit_min_enrichment: float = 1.0,
    bootstrap_replicates: int = 0,
    bootstrap_seed: int = 20260819,
    bootstrap_lower_quantile: float = 0.10,
    bootstrap_top_n: int = 50,
    bootstrap_min_probability_above_reference: float = 0.90,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    if pseudocount <= 0:
        raise ValueError("pseudocount must be positive")
    if not 0 < overall_gate_fraction <= 1:
        raise ValueError("overall_gate_fraction must be between 0 and 1")
    if strict_min_unsorted_count < 0 or strict_min_total_bin_count < 0:
        raise ValueError("strict count cutoffs cannot be negative")
    if not 0 <= strict_relative_median_fraction <= 1:
        raise ValueError("strict_relative_median_fraction must be between 0 and 1")
    if strict_min_high_bin_count < 0:
        raise ValueError("strict_min_high_bin_count cannot be negative")
    if not 1 <= strict_min_detected_bins <= 6:
        raise ValueError("strict_min_detected_bins must be between 1 and 6")
    if not 0 < jackpot_max_bin_probability <= 1:
        raise ValueError("jackpot_max_bin_probability must be between 0 and 1")
    if not 1 <= jackpot_max_detected_bins <= 6:
        raise ValueError("jackpot_max_detected_bins must be between 1 and 6")
    if min(
        top_hit_min_unsorted_count,
        top_hit_min_total_bin_count,
        top_hit_min_high_bin_count,
    ) < 0:
        raise ValueError("top-hit count cutoffs cannot be negative")
    if top_hit_min_enrichment <= 0:
        raise ValueError("top_hit_min_enrichment must be positive")
    if bootstrap_replicates < 0:
        raise ValueError("bootstrap_replicates cannot be negative")
    if not 0 < bootstrap_lower_quantile < 0.5:
        raise ValueError("bootstrap_lower_quantile must be between 0 and 0.5")
    if bootstrap_top_n <= 0:
        raise ValueError("bootstrap_top_n must be positive")
    if not 0 <= bootstrap_min_probability_above_reference <= 1:
        raise ValueError(
            "bootstrap_min_probability_above_reference must be between 0 and 1"
        )

    bins = samples[samples["sample_type"].eq("bin")].sort_values("bin_number")
    bin_ids = bins["sample_id"].astype(str).tolist()
    fractions = bins.set_index("sample_id").loc[bin_ids, "population_fraction"].astype(float)
    unsorted_id = str(samples.loc[samples["sample_type"].eq("unsorted"), "sample_id"].iloc[0])

    n_variants = len(variants)
    depths = counts.sum(axis=0)
    if (depths <= 0).any():
        raise ValueError(f"empty count samples: {depths[depths <= 0].index.tolist()}")

    raw_bins = counts[bin_ids].astype(float)
    raw_unsorted = counts[unsorted_id].astype(float)
    within_bin_frequency = raw_bins.div(depths[bin_ids], axis=1)
    unsorted_frequency = raw_unsorted / depths[unsorted_id]

    # This is the estimated amount of each UTR contributed by each population bin.
    population_mass = within_bin_frequency.mul(fractions, axis=1)
    reconstructed_gate_frequency = population_mass.sum(axis=1)
    bin_probability = population_mass.div(reconstructed_gate_frequency.replace(0, np.nan), axis=0)

    smooth_bins = (raw_bins + pseudocount).div(
        depths[bin_ids] + pseudocount * n_variants,
        axis=1,
    )
    smooth_unsorted = (raw_unsorted + pseudocount) / (
        depths[unsorted_id] + pseudocount * n_variants
    )
    smooth_gate = smooth_bins.mul(fractions, axis=1).sum(axis=1)

    total_bin_count = raw_bins.sum(axis=1)
    pass_coverage = (raw_unsorted >= min_unsorted_count) & (
        total_bin_count >= min_total_bin_count
    )
    detected_bin_count = (raw_bins > 0).sum(axis=1)

    result = variants.copy().set_index(id_column)
    result["unsorted_count"] = raw_unsorted
    result["total_6bin_count"] = total_bin_count
    result["detected_in_n_bins"] = detected_bin_count
    result["pass_coverage"] = pass_coverage
    result["unsorted_frequency"] = unsorted_frequency
    result["reconstructed_gate_frequency"] = reconstructed_gate_frequency

    bin_score = pd.Series(
        {sample_id: float(7 - int(number)) for sample_id, number in zip(bin_ids, bins["bin_number"])},
        dtype=float,
    )
    expected_bin_score = bin_probability.mul(bin_score, axis=1).sum(axis=1)
    baseline_score = float(np.dot(fractions.to_numpy(), bin_score.loc[bin_ids].to_numpy()))
    result["expected_bin_score"] = expected_bin_score
    result["neutral_baseline_score"] = baseline_score
    result["score_shift_from_pool"] = expected_bin_score - baseline_score

    for number, sample_id in enumerate(bin_ids, start=1):
        result[f"bin{number}_count"] = raw_bins[sample_id]
        result[f"bin{number}_probability"] = bin_probability[sample_id]
        result[f"bin{number}_vs_unsorted_log2"] = np.log2(
            smooth_bins[sample_id] / smooth_unsorted
        )

    high5_fraction = float(fractions.iloc[0])
    high15_fraction = float(fractions.iloc[:2].sum())
    low_fraction = float(fractions.iloc[4:].sum())
    result["high5_probability"] = bin_probability[bin_ids[0]]
    result["high15_probability"] = bin_probability[bin_ids[:2]].sum(axis=1)
    result["high15_enrichment"] = result["high15_probability"] / high15_fraction
    result["high15_log2_enrichment"] = np.log2(result["high15_enrichment"].clip(lower=1e-12))

    # Secondary/exploratory endpoint: combined high-bin composition relative to
    # whole unsorted. This includes both conditional High15 position and target-
    # gate representation, so it is retained for QC but not used as primary rank.
    smooth_high15_frequency = (
        smooth_bins[bin_ids[:2]].mul(fractions.iloc[:2], axis=1).sum(axis=1)
        / high15_fraction
    )
    bin1_vs_unsorted_enrichment = smooth_bins[bin_ids[0]] / smooth_unsorted
    bin2_vs_unsorted_enrichment = smooth_bins[bin_ids[1]] / smooth_unsorted
    top15_vs_unsorted_enrichment = smooth_high15_frequency / smooth_unsorted
    result["bin1_vs_unsorted_enrichment"] = bin1_vs_unsorted_enrichment
    result["bin2_vs_unsorted_enrichment"] = bin2_vs_unsorted_enrichment
    result["bin1_vs_unsorted_log2_enrichment"] = np.log2(
        bin1_vs_unsorted_enrichment.clip(lower=1e-12)
    )
    result["bin2_vs_unsorted_log2_enrichment"] = np.log2(
        bin2_vs_unsorted_enrichment.clip(lower=1e-12)
    )
    result["top15_combined_frequency"] = smooth_high15_frequency
    result["top15_vs_unsorted_enrichment"] = top15_vs_unsorted_enrichment
    result["top15_vs_unsorted_log2_enrichment"] = np.log2(
        top15_vs_unsorted_enrichment.clip(lower=1e-12)
    )
    result["low_bin_probability"] = bin_probability[bin_ids[4:]].sum(axis=1)
    result["low_bin_enrichment"] = result["low_bin_probability"] / low_fraction

    relative_bin_enrichment = bin_probability.div(fractions.to_numpy(), axis=1)
    enriched_number = relative_bin_enrichment.to_numpy().argmax(axis=1) + 1
    dominant_number = bin_probability.fillna(-1).to_numpy().argmax(axis=1) + 1
    result["most_enriched_bin"] = [f"bin{x}" for x in enriched_number]
    result["dominant_cell_mass_bin"] = [f"bin{x}" for x in dominant_number]
    probability_array = bin_probability.fillna(0).to_numpy()
    safe_probability_array = np.clip(probability_array, 1e-300, None)
    entropy = -np.sum(
        np.where(
            probability_array > 0,
            probability_array * np.log(safe_probability_array),
            0,
        ),
        axis=1,
    ) / math.log(6)
    result["distribution_entropy_0to1"] = entropy

    # Whole-population unsorted is used only for gate representation, not for the
    # within-gate fluorescence score.
    gate_representation_ratio = smooth_gate / smooth_unsorted
    gate_entry_raw = overall_gate_fraction * gate_representation_ratio
    result["gate_representation_ratio"] = gate_representation_ratio
    result["gate_entry_probability_raw"] = gate_entry_raw
    result["gate_entry_probability_capped"] = gate_entry_raw.clip(lower=0, upper=1)
    result["gate_entry_over_1_qc"] = gate_entry_raw > 1.10
    high15_population_mass = population_mass[bin_ids[:2]].sum(axis=1)
    result["whole_population_high15_probability_raw"] = (
        overall_gate_fraction * high15_population_mass / smooth_unsorted
    )

    mfi = bins.set_index("sample_id").loc[bin_ids, "representative_log10_mfi"]
    if mfi.notna().all():
        result["expected_log10_mfi"] = bin_probability.mul(mfi.astype(float), axis=1).sum(axis=1)

    plasmid_rows = samples[samples["sample_type"].eq("plasmid")]
    if len(plasmid_rows) == 1:
        plasmid_id = str(plasmid_rows["sample_id"].iloc[0])
        raw_plasmid = counts[plasmid_id].astype(float)
        smooth_plasmid = (raw_plasmid + pseudocount) / (
            depths[plasmid_id] + pseudocount * n_variants
        )
        result["plasmid_count"] = raw_plasmid
        result["unsorted_vs_plasmid_log2"] = np.log2(smooth_unsorted / smooth_plasmid)

    score_columns = [
        "expected_bin_score",
        "score_shift_from_pool",
        "high5_probability",
        "high15_probability",
        "high15_enrichment",
        "high15_log2_enrichment",
        "bin1_vs_unsorted_enrichment",
        "bin2_vs_unsorted_enrichment",
        "bin1_vs_unsorted_log2_enrichment",
        "bin2_vs_unsorted_log2_enrichment",
        "top15_combined_frequency",
        "top15_vs_unsorted_enrichment",
        "top15_vs_unsorted_log2_enrichment",
        "low_bin_probability",
        "low_bin_enrichment",
    ]
    if "expected_log10_mfi" in result.columns:
        score_columns.append("expected_log10_mfi")
    result.loc[~pass_coverage, score_columns] = np.nan
    result.loc[~pass_coverage, ["most_enriched_bin", "dominant_cell_mass_bin"]] = "low_coverage"

    result["estimated_rank"] = result["expected_bin_score"].rank(
        ascending=False, method="average"
    )
    passing_n = int(result["expected_bin_score"].notna().sum())
    if passing_n > 1:
        result["expression_percentile"] = 100 * (
            1 - (result["estimated_rank"] - 1) / (passing_n - 1)
        )
    elif passing_n == 1:
        result["expression_percentile"] = np.where(
            result["expected_bin_score"].notna(), 100.0, np.nan
        )
    else:
        result["expression_percentile"] = np.nan
    result["expression_tier"] = result["expression_percentile"].map(percentile_tier)
    result["high_candidate_flag"] = (
        result["pass_coverage"]
        & (result["expression_percentile"] >= 95)
        & (result["high15_enrichment"] >= 1.5)
    )

    # Strict candidate-confidence filter. The fixed floors prevent sparse reads
    # from dominating the score, while the median-relative term scales with the
    # sequencing depth of each experiment. This is intentionally separate from
    # the permissive pass_coverage filter used to preserve the full exploration
    # table.
    median_unsorted_count = float(raw_unsorted.median())
    median_total_bin_count = float(total_bin_count.median())
    strict_unsorted_cutoff = max(
        int(strict_min_unsorted_count),
        int(round(strict_relative_median_fraction * median_unsorted_count)),
    )
    strict_total_bin_cutoff = max(
        int(strict_min_total_bin_count),
        int(round(strict_relative_median_fraction * median_total_bin_count)),
    )
    result["high_bin_raw_count"] = raw_bins.iloc[:, :2].sum(axis=1)
    result["maximum_bin_probability"] = bin_probability.max(axis=1)
    result["strict_unsorted_cutoff"] = strict_unsorted_cutoff
    result["strict_total_6bin_cutoff"] = strict_total_bin_cutoff
    result["strict_coverage_pass"] = (
        (result["unsorted_count"] >= strict_unsorted_cutoff)
        & (result["total_6bin_count"] >= strict_total_bin_cutoff)
        & (result["high_bin_raw_count"] >= strict_min_high_bin_count)
        & (result["detected_in_n_bins"] >= strict_min_detected_bins)
    )
    result["single_bin_jackpot_suspect"] = (
        (result["maximum_bin_probability"] >= jackpot_max_bin_probability)
        & (result["detected_in_n_bins"] <= jackpot_max_detected_bins)
    )

    # Primary endpoint: conditional probability of occupying bin1 or bin2 among
    # cells already inside the mCherry+/GFP- target gate.  Unsorted is retained
    # for read support and gate-representation QC, but it does not define this
    # fluorescence ranking.  A real sharp shift can concentrate in bin1 and
    # deplete bin2, so individual enrichment in both bins is diagnostic only.
    result["top15_read_support_pass"] = (
        (result["unsorted_count"] >= top_hit_min_unsorted_count)
        & (result["total_6bin_count"] >= top_hit_min_total_bin_count)
        & (result["high_bin_raw_count"] >= top_hit_min_high_bin_count)
        & result["high15_probability"].notna()
    )
    result["top15_both_bins_enriched"] = (
        (result["bin1_vs_unsorted_enrichment"] >= top_hit_min_enrichment)
        & (result["bin2_vs_unsorted_enrichment"] >= top_hit_min_enrichment)
    )
    high15_rank_values = result["high15_probability"].where(
        result["top15_read_support_pass"]
    )
    result["high15_primary_rank"] = high15_rank_values.rank(
        ascending=False, method="average"
    )
    high15_rank_n = int(high15_rank_values.notna().sum())
    if high15_rank_n > 1:
        result["high15_primary_percentile"] = 100 * (
            1 - (result["high15_primary_rank"] - 1) / (high15_rank_n - 1)
        )
    elif high15_rank_n == 1:
        result["high15_primary_percentile"] = np.where(
            high15_rank_values.notna(), 100.0, np.nan
        )
    else:
        result["high15_primary_percentile"] = np.nan
    high15_comparator = high15_fraction * top_hit_min_enrichment
    result["high15_comparator_probability"] = high15_comparator
    result["high15_log2_fold_vs_comparator"] = np.log2(
        result["high15_probability"].clip(lower=1e-12) / high15_comparator
    )
    result["top15_above_comparator"] = result["high15_probability"] > high15_comparator
    result["high_confidence_candidate_flag"] = (
        result["strict_coverage_pass"]
        & result["high_candidate_flag"]
        & ~result["single_bin_jackpot_suspect"]
    )

    strict_scores = result["expected_bin_score"].where(result["strict_coverage_pass"])
    result["strict_estimated_rank"] = strict_scores.rank(ascending=False, method="average")
    strict_n = int(strict_scores.notna().sum())
    if strict_n > 1:
        result["strict_expression_percentile"] = 100 * (
            1 - (result["strict_estimated_rank"] - 1) / (strict_n - 1)
        )
    elif strict_n == 1:
        result["strict_expression_percentile"] = np.where(
            strict_scores.notna(), 100.0, np.nan
        )
    else:
        result["strict_expression_percentile"] = np.nan
    result["strict_expression_tier"] = result["strict_expression_percentile"].map(
        percentile_tier
    )

    reference_summary: dict[str, object] = {
        "reference_variant_id": reference_variant_id,
        "reference_display_label": reference_label,
        "reference_comparison_enabled": reference_variant_id is not None,
    }
    reference_score = baseline_score
    reference_high15 = high15_comparator
    result["is_reference_variant"] = False
    if reference_variant_id is not None:
        if reference_variant_id not in result.index:
            raise ValueError(f"reference variant ID not found in analysis result: {reference_variant_id!r}")
        reference = result.loc[reference_variant_id]
        if not bool(reference["pass_coverage"]):
            raise ValueError(
                f"reference variant {reference_variant_id!r} failed the coverage filter; "
                "reference-relative ranking would be unreliable"
            )
        reference_score = float(reference["expected_bin_score"])
        reference_high15 = float(reference["high15_probability"])
        reference_high15_enrichment = float(reference["high15_enrichment"])
        reference_top15_unsorted_enrichment = float(
            reference["top15_vs_unsorted_enrichment"]
        )
        reference_top15_unsorted_log2 = float(
            reference["top15_vs_unsorted_log2_enrichment"]
        )
        result.loc[reference_variant_id, "is_reference_variant"] = True
        result["reference_display_label"] = reference_label or reference_variant_id
        result["reference_expected_bin_score"] = reference_score
        result["delta_score_vs_reference"] = result["expected_bin_score"] - reference_score
        result["reference_high15_probability"] = reference_high15
        result["delta_high15_probability_vs_reference"] = (
            result["high15_probability"] - reference_high15
        )
        result["high15_fold_vs_reference"] = (
            result["high15_probability"] / reference_high15
            if reference_high15 > 0
            else np.nan
        )
        result["high15_log2_fold_vs_reference"] = np.log2(
            result["high15_fold_vs_reference"].clip(lower=1e-12)
        )
        result["high15_comparator_probability"] = reference_high15
        result["high15_log2_fold_vs_comparator"] = result[
            "high15_log2_fold_vs_reference"
        ]
        result["score_above_reference"] = (
            result["pass_coverage"] & (result["delta_score_vs_reference"] > 0)
        )
        result["score_and_high15_above_reference"] = (
            result["score_above_reference"]
            & (result["delta_high15_probability_vs_reference"] > 0)
        )
        result["reference_top15_vs_unsorted_enrichment"] = (
            reference_top15_unsorted_enrichment
        )
        result["reference_top15_vs_unsorted_log2_enrichment"] = (
            reference_top15_unsorted_log2
        )
        result["delta_top15_log2_enrichment_vs_reference"] = (
            result["top15_vs_unsorted_log2_enrichment"]
            - reference_top15_unsorted_log2
        )
        result["top15_enrichment_fold_vs_reference"] = (
            result["top15_vs_unsorted_enrichment"]
            / reference_top15_unsorted_enrichment
            if reference_top15_unsorted_enrichment > 0
            else np.nan
        )
        result["top15_above_comparator"] = (
            result["delta_high15_probability_vs_reference"] > 0
        )
        reference_summary.update(
            {
                "reference_pass_coverage": True,
                "reference_expected_bin_score": reference_score,
                "reference_high15_probability": reference_high15,
                "reference_high15_enrichment": reference_high15_enrichment,
                "reference_top15_vs_unsorted_enrichment": reference_top15_unsorted_enrichment,
                "reference_top15_vs_unsorted_log2_enrichment": reference_top15_unsorted_log2,
                "reference_estimated_rank": float(reference["estimated_rank"]),
                "reference_expression_percentile": float(reference["expression_percentile"]),
                "variants_with_score_above_reference": int(
                    result["score_above_reference"].sum()
                ),
                "variants_with_score_and_high15_above_reference": int(
                    result["score_and_high15_above_reference"].sum()
                ),
                "high_candidates_above_reference": int(
                    (
                        result["high_candidate_flag"]
                        & result["score_and_high15_above_reference"]
                    ).sum()
                ),
            }
        )

    bootstrap = technical_high15_bootstrap(
        raw_bins=raw_bins,
        fractions=fractions,
        reference_variant_id=reference_variant_id,
        read_support=result["top15_read_support_pass"],
        replicates=bootstrap_replicates,
        seed=bootstrap_seed,
        lower_quantile=bootstrap_lower_quantile,
        top_n=bootstrap_top_n,
    )
    result = result.join(bootstrap)
    result["high15_robust_rank_score"] = result[
        "high15_log2_fold_bootstrap_lower"
    ]
    if bootstrap_replicates > 0:
        final_rank_values = result["high15_robust_rank_score"].where(
            result["top15_read_support_pass"]
        )
        result["high15_technical_stability_pass"] = (
            result["high15_bootstrap_probability_above_comparator"]
            >= bootstrap_min_probability_above_reference
        )
    else:
        final_rank_values = result["high15_probability"].where(
            result["top15_read_support_pass"]
        )
        result["high15_technical_stability_pass"] = result[
            "top15_above_comparator"
        ]
    result["high15_final_rank"] = final_rank_values.rank(
        ascending=False, method="average"
    )
    final_rank_n = int(final_rank_values.notna().sum())
    if final_rank_n > 1:
        result["high15_final_percentile"] = 100 * (
            1 - (result["high15_final_rank"] - 1) / (final_rank_n - 1)
        )
    elif final_rank_n == 1:
        result["high15_final_percentile"] = np.where(
            final_rank_values.notna(), 100.0, np.nan
        )
    else:
        result["high15_final_percentile"] = np.nan

    # Backward-compatible names now point to the conditional High15 primary
    # ranking.  The old unsorted-based endpoint remains available explicitly as
    # top15_vs_unsorted_* and is secondary/exploratory.
    result["top15_rank"] = result["high15_final_rank"]
    result["top15_percentile"] = result["high15_final_percentile"]
    result["high15_weighted_score_support"] = (
        result["expected_bin_score"] > reference_score
    )
    result["top15_candidate_flag"] = (
        result["top15_read_support_pass"] & result["top15_above_comparator"]
    )
    result["top15_priority_candidate_flag"] = (
        result["top15_candidate_flag"]
        & result["high15_technical_stability_pass"]
        & ~result["single_bin_jackpot_suspect"]
    )
    result["candidate_tier"] = np.select(
        [
            result["top15_priority_candidate_flag"]
            & result["high15_weighted_score_support"],
            result["top15_priority_candidate_flag"],
            result["top15_candidate_flag"],
        ],
        ["tier1_clean_high_shift", "tier2_high_tail", "tier3_review"],
        default="not_candidate",
    )
    result["recommended_for_cloning"] = result["candidate_tier"].isin(
        {"tier1_clean_high_shift", "tier2_high_tail"}
    )
    result["default_top_n_cloning_shortlist"] = (
        result["recommended_for_cloning"]
        & (result["high15_final_rank"] <= bootstrap_top_n)
        & ~result["is_reference_variant"]
    )

    result = result.reset_index().sort_values(
        ["top15_read_support_pass", "high15_final_rank", "estimated_rank"],
        ascending=[False, True, True],
        na_position="last",
    )

    essential_columns = [
        id_column,
        "pass_coverage",
        "unsorted_count",
        "total_6bin_count",
        "expected_bin_score",
        "score_shift_from_pool",
        "high5_probability",
        "high15_probability",
        "high15_enrichment",
        "high15_comparator_probability",
        "high15_log2_fold_vs_comparator",
        "high15_primary_rank",
        "high15_primary_percentile",
        "high15_robust_rank_score",
        "high15_final_rank",
        "high15_final_percentile",
        "high15_bootstrap_probability_above_comparator",
        "high15_bootstrap_top_n_frequency",
        "high15_technical_stability_pass",
        "high15_weighted_score_support",
        "bin1_vs_unsorted_enrichment",
        "bin2_vs_unsorted_enrichment",
        "top15_vs_unsorted_enrichment",
        "top15_vs_unsorted_log2_enrichment",
        "top15_read_support_pass",
        "top15_both_bins_enriched",
        "top15_rank",
        "top15_percentile",
        "top15_candidate_flag",
        "top15_priority_candidate_flag",
        "candidate_tier",
        "recommended_for_cloning",
        "default_top_n_cloning_shortlist",
        "most_enriched_bin",
        "gate_entry_probability_capped",
        "estimated_rank",
        "expression_percentile",
        "expression_tier",
        "high_candidate_flag",
        "strict_coverage_pass",
        "high_bin_raw_count",
        "maximum_bin_probability",
        "single_bin_jackpot_suspect",
        "strict_estimated_rank",
        "strict_expression_percentile",
        "strict_expression_tier",
        "high_confidence_candidate_flag",
    ]
    if reference_variant_id is not None:
        essential_columns.extend(
            [
                "is_reference_variant",
                "delta_score_vs_reference",
                "delta_high15_probability_vs_reference",
                "high15_fold_vs_reference",
                "high15_log2_fold_vs_reference",
                "score_above_reference",
                "score_and_high15_above_reference",
                "delta_top15_log2_enrichment_vs_reference",
                "top15_enrichment_fold_vs_reference",
            ]
        )
    if "expected_log10_mfi" in result.columns:
        essential_columns.insert(7, "expected_log10_mfi")
    essential = result[essential_columns].copy()

    valid = result[result["pass_coverage"]].copy()
    log_unsorted = np.log10(valid["unsorted_frequency"].clip(lower=1e-15))
    log_gate = np.log10(valid["reconstructed_gate_frequency"].clip(lower=1e-15))
    if len(valid) >= 3 and log_unsorted.nunique() > 1 and log_gate.nunique() > 1:
        correlation = spearmanr(
            log_unsorted,
            log_gate,
        )
        gate_unsorted_rho = float(correlation.statistic)
        gate_unsorted_p = float(correlation.pvalue)
    else:
        gate_unsorted_rho = float("nan")
        gate_unsorted_p = float("nan")
    fraction_source = (
        str(bins["population_fraction_source"].iloc[0])
        if "population_fraction_source" in bins.columns
        else "population_fraction"
    )
    summary = {
        "reference_variants": int(len(variants)),
        "variants_passing_coverage": int(pass_coverage.sum()),
        "high_candidate_count": int(result["high_candidate_flag"].sum()),
        "strict_coverage_passing_count": int(result["strict_coverage_pass"].sum()),
        "single_bin_jackpot_suspect_count": int(
            result["single_bin_jackpot_suspect"].sum()
        ),
        "high_confidence_candidate_count": int(
            result["high_confidence_candidate_flag"].sum()
        ),
        "top15_read_support_passing_count": int(
            result["top15_read_support_pass"].sum()
        ),
        "top15_candidate_count": int(result["top15_candidate_flag"].sum()),
        "top15_priority_candidate_count": int(
            result["top15_priority_candidate_flag"].sum()
        ),
        "tier1_clean_high_shift_count": int(
            result["candidate_tier"].eq("tier1_clean_high_shift").sum()
        ),
        "tier2_high_tail_count": int(
            result["candidate_tier"].eq("tier2_high_tail").sum()
        ),
        "tier3_review_count": int(
            result["candidate_tier"].eq("tier3_review").sum()
        ),
        "recommended_for_cloning_count": int(
            result["recommended_for_cloning"].sum()
        ),
        "default_top_n_cloning_shortlist_count": int(
            result["default_top_n_cloning_shortlist"].sum()
        ),
        "primary_endpoint": "conditional_high15_probability",
        "secondary_endpoint": "expected_bin_score_6_to_1",
        "unsorted_endpoint_role": "gate_representation_qc_and_exploratory",
        "top_hit_min_unsorted_count": top_hit_min_unsorted_count,
        "top_hit_min_total_bin_count": top_hit_min_total_bin_count,
        "top_hit_min_high_bin_count": top_hit_min_high_bin_count,
        "top_hit_min_enrichment": top_hit_min_enrichment,
        "bootstrap_replicates": bootstrap_replicates,
        "bootstrap_seed": bootstrap_seed,
        "bootstrap_lower_quantile": bootstrap_lower_quantile,
        "bootstrap_top_n": bootstrap_top_n,
        "bootstrap_min_probability_above_reference": (
            bootstrap_min_probability_above_reference
        ),
        "bootstrap_uncertainty_scope": "technical_ngs_read_sampling_only",
        "median_unsorted_count": median_unsorted_count,
        "median_total_6bin_count": median_total_bin_count,
        "strict_unsorted_cutoff": strict_unsorted_cutoff,
        "strict_total_6bin_cutoff": strict_total_bin_cutoff,
        "strict_min_high_bin_count": strict_min_high_bin_count,
        "strict_min_detected_bins": strict_min_detected_bins,
        "jackpot_max_bin_probability": jackpot_max_bin_probability,
        "jackpot_max_detected_bins": jackpot_max_detected_bins,
        "neutral_baseline_score": baseline_score,
        "population_fraction_source": fraction_source,
        "overall_gate_fraction": overall_gate_fraction,
        "gate_vs_unsorted_spearman_rho": gate_unsorted_rho,
        "gate_vs_unsorted_spearman_p": gate_unsorted_p,
        "min_unsorted_count": min_unsorted_count,
        "min_total_bin_count": min_total_bin_count,
        **reference_summary,
    }
    return result, essential, summary


def make_plots(
    outdir: Path,
    sample_qc: pd.DataFrame,
    result: pd.DataFrame,
    id_column: str,
    top_n: int,
) -> list[Path]:
    sns.set_theme(style="whitegrid")
    plot_dir = outdir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    fig, ax = plt.subplots(figsize=(8, 4.5))
    sns.barplot(
        data=sample_qc,
        x="sample_id",
        y="assigned_reference_reads",
        hue="sample_type",
        dodge=False,
        ax=ax,
    )
    ax.set_yscale("log")
    ax.set_xlabel("")
    ax.set_ylabel("Assigned reads (log scale)")
    ax.tick_params(axis="x", rotation=35)
    fig.tight_layout()
    path = plot_dir / "01_sample_read_depth.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    written.append(path)

    passing = result[result["pass_coverage"]].copy()
    reference_rows = (
        passing[passing["is_reference_variant"]]
        if "is_reference_variant" in passing.columns
        else passing.iloc[0:0]
    )
    if not passing.empty:
        fig, ax = plt.subplots(figsize=(6.2, 4.8))
        sns.scatterplot(
            data=passing,
            x="expected_bin_score",
            y="high15_enrichment",
            hue="expression_tier",
            size="unsorted_count",
            sizes=(12, 90),
            alpha=0.7,
            ax=ax,
        )
        ax.axhline(1.0, color="black", linestyle="--", linewidth=1)
        if len(reference_rows) == 1:
            reference = reference_rows.iloc[0]
            ax.scatter(
                reference["expected_bin_score"],
                reference["high15_enrichment"],
                marker="*",
                s=260,
                color="#D62728",
                edgecolor="white",
                linewidth=0.9,
                zorder=10,
            )
            ax.annotate(
                f"reference: {reference.get('reference_display_label', reference[id_column])}",
                (reference["expected_bin_score"], reference["high15_enrichment"]),
                xytext=(8, 8),
                textcoords="offset points",
                color="#B2182B",
                weight="bold",
            )
        ax.set_xlabel("Expected bin score (6=high, 1=low)")
        ax.set_ylabel("Top-15% enrichment (1=pool average)")
        fig.tight_layout()
        path = plot_dir / "02_expression_score_vs_high15.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        written.append(path)

        top = passing.nsmallest(min(top_n, len(passing)), "estimated_rank")
        if len(reference_rows) == 1 and not top[id_column].eq(reference_rows.iloc[0][id_column]).any():
            top = pd.concat([top, reference_rows], ignore_index=True)
        probability_columns = [f"bin{x}_probability" for x in range(1, 7)]
        display_ids = top[id_column].astype(str).where(
            ~top["is_reference_variant"],
            top.get("reference_display_label", top[id_column]).astype(str) + " [reference]",
        )
        profile = top.assign(_display_id=display_ids).set_index("_display_id")[probability_columns]
        fig, ax = plt.subplots(figsize=(7.2, max(4.5, min(16, len(profile) * 0.24))))
        sns.heatmap(
            profile,
            cmap="mako",
            vmin=0,
            vmax=max(0.25, float(profile.max().max())),
            cbar_kws={"label": "Estimated fraction of this UTR's gated cells"},
            ax=ax,
        )
        ax.set_xlabel("bin1 = highest mCherry; bin6 = lowest mCherry-positive")
        ax.set_ylabel(id_column)
        ax.set_title(
            f"Top {top_n} UTR bin distributions"
            + (" + reference" if len(reference_rows) == 1 else "")
        )
        fig.tight_layout()
        path = plot_dir / "03_top_utr_bin_heatmap.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        written.append(path)

        positive = passing[
            (passing["unsorted_frequency"] > 0)
            & (passing["reconstructed_gate_frequency"] > 0)
        ]
        if not positive.empty:
            fig, ax = plt.subplots(figsize=(5.3, 5.0))
            ax.scatter(
                np.log10(positive["unsorted_frequency"]),
                np.log10(positive["reconstructed_gate_frequency"]),
                s=10,
                alpha=0.45,
            )
            reference_positive = positive[positive["is_reference_variant"]]
            if len(reference_positive) == 1:
                reference = reference_positive.iloc[0]
                ref_x = np.log10(reference["unsorted_frequency"])
                ref_y = np.log10(reference["reconstructed_gate_frequency"])
                ax.scatter(
                    ref_x,
                    ref_y,
                    marker="*",
                    s=260,
                    color="#D62728",
                    edgecolor="white",
                    linewidth=0.9,
                    zorder=10,
                )
                ax.annotate(
                    f"reference: {reference.get('reference_display_label', reference[id_column])}",
                    (ref_x, ref_y),
                    xytext=(8, 8),
                    textcoords="offset points",
                    color="#B2182B",
                    weight="bold",
                )
            lo = min(ax.get_xlim()[0], ax.get_ylim()[0])
            hi = max(ax.get_xlim()[1], ax.get_ylim()[1])
            ax.plot([lo, hi], [lo, hi], "--", color="black", linewidth=1)
            ax.set_xlabel("log10 UTR frequency in whole unsorted")
            ax.set_ylabel("log10 reconstructed frequency in target gate")
            ax.set_title("Whole unsorted vs mCherry+/GFP- gate")
            fig.tight_layout()
            path = plot_dir / "04_unsorted_vs_target_gate.png"
            fig.savefig(path, dpi=180)
            plt.close(fig)
            written.append(path)

        if len(reference_rows) == 1:
            reference = reference_rows.iloc[0]
            fig, ax = plt.subplots(figsize=(6.2, 4.5))
            ax.hist(passing["expected_bin_score"], bins=40, color="#4C78A8", alpha=0.8)
            ax.axvline(
                reference["expected_bin_score"],
                color="#D62728",
                linewidth=2.2,
                linestyle="--",
                label=(
                    f"{reference.get('reference_display_label', reference[id_column])}: "
                    f"{reference['expected_bin_score']:.3f} "
                    f"(rank {reference['estimated_rank']:.0f})"
                ),
            )
            ax.set_xlabel("Expected bin score")
            ax.set_ylabel("Coverage-passing UTRs")
            ax.set_title("Reference position in the fluorescence-score distribution")
            ax.legend(frameon=False)
            fig.tight_layout()
            path = plot_dir / "05_reference_score_position.png"
            fig.savefig(path, dpi=180)
            plt.close(fig)
            written.append(path)
    return written


def dataframe_html(table: pd.DataFrame, rows: int = 30) -> str:
    shown = table.head(rows).copy()
    for column in shown.select_dtypes(include="number").columns:
        shown[column] = shown[column].map(
            lambda value: "" if pd.isna(value) else f"{value:.4g}"
        )
    return shown.to_html(index=False, escape=True, border=0)


def write_report(
    outdir: Path,
    sample_qc: pd.DataFrame,
    essential: pd.DataFrame,
    summary: dict[str, object],
    plots: list[Path],
    id_column: str,
) -> Path:
    report = outdir / "easy_report.html"
    images = "\n".join(
        f'<figure><img src="{html.escape(str(path.relative_to(outdir)))}" '
        f'alt="{html.escape(path.stem)}"><figcaption>{html.escape(path.stem)}</figcaption></figure>'
        for path in plots
    )
    if summary.get("reference_comparison_enabled"):
        reference_html = (
            '<div class="box"><b>Reference 위치</b>: '
            f"{html.escape(str(summary.get('reference_display_label') or summary['reference_variant_id']))}; "
            f"score {float(summary['reference_expected_bin_score']):.3f}; "
            f"rank {float(summary['reference_estimated_rank']):.0f}; "
            f"percentile {float(summary['reference_expression_percentile']):.2f}. "
            f"Reference보다 score가 높은 UTR은 "
            f"{int(summary['variants_with_score_above_reference']):,}개입니다.</div>"
        )
    else:
        reference_html = ""
    report.write_text(
        f"""<!doctype html><html lang="ko"><head><meta charset="utf-8">
<title>7-sample 5' UTR Sort-seq report</title><style>
body{{font-family:Arial,'Noto Sans KR',sans-serif;max-width:1200px;margin:26px auto;padding:0 20px;color:#222}}
h1,h2{{color:#17365d}} .box{{background:#eef5fb;border-left:5px solid #4472c4;padding:12px 16px}}
.warn{{background:#fff2cc;border-left-color:#bf9000}} table{{border-collapse:collapse;font-size:12px;display:block;overflow:auto}}
th,td{{border:1px solid #ddd;padding:5px 7px}} th{{background:#d9eaf7}} img{{max-width:100%;height:auto}}
figure{{margin:20px 0}} code{{background:#f2f2f2;padding:2px 4px}}
</style></head><body><h1>7-sample 5' UTR Sort-seq 분석</h1>
<div class="box"><b>한 줄 해석</b>: 각 UTR의 reads를 sample depth로 나눈 뒤 bin 크기를 곱해,
그 UTR 세포가 bin1–6에 어떻게 분포하는지 복원했습니다. bin1은 최고 mCherry, bin6은 최저 mCherry-positive입니다.</div>
<p>Reference UTR: {summary['reference_variants']:,}; coverage 통과: {summary['variants_passing_coverage']:,};
high-candidate heuristic 통과: {summary['high_candidate_count']:,};
strict coverage 통과: {summary['strict_coverage_passing_count']:,};
high-confidence 후보: {summary['high_confidence_candidate_count']:,}</p>
{reference_html}
<h2>결과를 읽는 순서</h2><ol>
<li><code>pass_coverage</code>가 TRUE인 UTR만 봅니다.</li>
<li><code>high15_final_rank</code>가 cloning 후보의 1차 순위입니다. 낮을수록 우선입니다.</li>
<li><code>high15_probability</code>는 target gate 안의 해당 UTR 세포 중 bin1+2에 있을 추정 비율입니다.</li>
<li><code>expected_bin_score</code>는 전체 6-bin 이동을 확인하는 보조 지표입니다.</li>
<li><code>bin1_probability</code>–<code>bin6_probability</code>로 분포가 자연스러운지 확인합니다.</li>
<li><code>candidate_tier</code>, 기술적 bootstrap 안정성, jackpot 경고를 함께 확인합니다.</li>
</ol>
<div class="box warn"><b>주의</b>: biological replicate가 없는 한 이것은 후보 선별용 추정 순위입니다.
bootstrap은 NGS read sampling만 반영하며 PCR/생물학적 불확실성을 대체하지 않습니다. FDR 유의성이나 translation rate의 직접 측정값으로 해석하지 마세요.</div>
<h2>Sample QC</h2>{dataframe_html(sample_qc, 20)}
<h2>상위 UTR 요약</h2>{dataframe_html(essential, 50)}
<h2>Plots</h2>{images}
</body></html>""",
        encoding="utf-8",
    )
    return report


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    variants = load_variants(args.variants, args.id_column)
    reference_variant_id = resolve_reference_variant_id(
        variants,
        args.id_column,
        args.reference_variant_id,
    )
    reference_label = reference_display_label(
        variants,
        args.id_column,
        reference_variant_id,
    )
    if args.reference_variant_id.lower() == "auto" and reference_variant_id is None:
        print(
            "WARNING: no original/orginal reference variant was found; "
            "reference comparison is disabled",
            file=sys.stderr,
        )
    samples = load_samples(args.samples)
    counts, sample_qc = load_count_matrix(samples, variants, args.id_column, args.samples)
    result, essential, summary = analyze(
        variants,
        samples,
        counts,
        args.id_column,
        args.pseudocount,
        args.min_unsorted_count,
        args.min_total_bin_count,
        args.overall_gate_fraction,
        reference_variant_id,
        reference_label,
        args.strict_min_unsorted_count,
        args.strict_min_total_bin_count,
        args.strict_relative_median_fraction,
        args.strict_min_high_bin_count,
        args.strict_min_detected_bins,
        args.jackpot_max_bin_probability,
        args.jackpot_max_detected_bins,
        args.top_hit_min_unsorted_count,
        args.top_hit_min_total_bin_count,
        args.top_hit_min_high_bin_count,
        args.top_hit_min_enrichment,
        args.bootstrap_replicates,
        args.bootstrap_seed,
        args.bootstrap_lower_quantile,
        args.bootstrap_top_n,
        args.bootstrap_min_probability_above_reference,
    )
    plots = make_plots(args.outdir, sample_qc, result, args.id_column, args.top_n)
    report = write_report(args.outdir, sample_qc, essential, summary, plots, args.id_column)

    result.to_csv(args.outdir / "utr_results_full.tsv", sep="\t", index=False)
    essential.to_csv(args.outdir / "utr_results_easy.tsv", sep="\t", index=False)
    essential[essential["high_candidate_flag"]].to_csv(
        args.outdir / "high_candidates.tsv", sep="\t", index=False
    )
    result[result["strict_coverage_pass"]].to_csv(
        args.outdir / "strict_coverage_results.tsv", sep="\t", index=False
    )
    result[result["high_confidence_candidate_flag"]].to_csv(
        args.outdir / "high_confidence_candidates.tsv", sep="\t", index=False
    )
    top15_export_columns = [
        args.id_column,
        "high15_final_rank",
        "high15_final_percentile",
        "high15_primary_rank",
        "high15_probability",
        "high15_enrichment",
        "high15_comparator_probability",
        "high15_log2_fold_vs_comparator",
        "high15_robust_rank_score",
        "high15_bootstrap_median",
        "high15_bootstrap_lower",
        "high15_bootstrap_upper",
        "high15_log2_fold_bootstrap_median",
        "high15_log2_fold_bootstrap_lower",
        "high15_log2_fold_bootstrap_upper",
        "high15_bootstrap_probability_above_comparator",
        "high15_bootstrap_top_n_frequency",
        "high15_technical_stability_pass",
        "unsorted_count",
        "total_6bin_count",
        "high_bin_raw_count",
        "detected_in_n_bins",
        "expected_bin_score",
        "high15_weighted_score_support",
        *[f"bin{x}_count" for x in range(1, 7)],
        *[f"bin{x}_probability" for x in range(1, 7)],
        "gate_representation_ratio",
        "top15_vs_unsorted_enrichment",
        "top15_vs_unsorted_log2_enrichment",
        "bin1_vs_unsorted_enrichment",
        "bin2_vs_unsorted_enrichment",
        "top15_both_bins_enriched",
        "single_bin_jackpot_suspect",
        "top15_candidate_flag",
        "top15_priority_candidate_flag",
        "candidate_tier",
        "recommended_for_cloning",
        "default_top_n_cloning_shortlist",
        "is_reference_variant",
    ]
    if reference_variant_id is not None:
        top15_export_columns.extend(
            [
                "reference_display_label",
                "delta_high15_probability_vs_reference",
                "high15_fold_vs_reference",
                "high15_log2_fold_vs_reference",
                "delta_top15_log2_enrichment_vs_reference",
                "top15_enrichment_fold_vs_reference",
            ]
        )
    top15_ranking = result.loc[
        result["top15_read_support_pass"], top15_export_columns
    ].sort_values("high15_final_rank", ascending=True, na_position="last")
    top15_ranking.to_csv(
        args.outdir / "high15_primary_ranking.csv",
        index=False,
        encoding="utf-8-sig",
    )
    # Backward-compatible filename retained for existing offline workflows.
    top15_ranking.to_csv(
        args.outdir / "top15_enrichment_ranking.csv",
        index=False,
        encoding="utf-8-sig",
    )
    top15_ranking[top15_ranking["top15_candidate_flag"]].to_csv(
        args.outdir / "high15_candidates.csv",
        index=False,
        encoding="utf-8-sig",
    )
    top15_ranking[top15_ranking["top15_candidate_flag"]].to_csv(
        args.outdir / "top15_candidates.csv",
        index=False,
        encoding="utf-8-sig",
    )
    top15_ranking[top15_ranking["top15_priority_candidate_flag"]].to_csv(
        args.outdir / "top15_priority_candidates.csv",
        index=False,
        encoding="utf-8-sig",
    )
    top15_ranking[
        top15_ranking["recommended_for_cloning"]
        & ~top15_ranking["is_reference_variant"]
    ].to_csv(
        args.outdir / "top_candidates_for_cloning.csv",
        index=False,
        encoding="utf-8-sig",
    )
    top15_ranking[top15_ranking["default_top_n_cloning_shortlist"]].to_csv(
        args.outdir / f"top{args.bootstrap_top_n}_candidates_for_cloning.csv",
        index=False,
        encoding="utf-8-sig",
    )
    if reference_variant_id is not None:
        comparison_columns = [
            args.id_column,
            "pass_coverage",
            "unsorted_count",
            "total_6bin_count",
            "top15_read_support_pass",
            "expected_bin_score",
            "delta_score_vs_reference",
            "high15_probability",
            "high15_final_rank",
            "high15_robust_rank_score",
            "high15_bootstrap_probability_above_comparator",
            "high15_bootstrap_top_n_frequency",
            "high15_technical_stability_pass",
            "high15_weighted_score_support",
            "candidate_tier",
            "recommended_for_cloning",
            "default_top_n_cloning_shortlist",
            "top15_vs_unsorted_enrichment",
            "top15_vs_unsorted_log2_enrichment",
            "delta_top15_log2_enrichment_vs_reference",
            "top15_enrichment_fold_vs_reference",
            "bin1_vs_unsorted_enrichment",
            "bin2_vs_unsorted_enrichment",
            "top15_rank",
            "top15_candidate_flag",
            "top15_priority_candidate_flag",
            "delta_high15_probability_vs_reference",
            "high15_fold_vs_reference",
            "high15_log2_fold_vs_reference",
            "estimated_rank",
            "expression_percentile",
            "expression_tier",
            "high_candidate_flag",
            "strict_coverage_pass",
            "single_bin_jackpot_suspect",
            "strict_estimated_rank",
            "strict_expression_percentile",
            "high_confidence_candidate_flag",
            "is_reference_variant",
            "reference_display_label",
        ]
        result.loc[:, comparison_columns].sort_values(
            ["top15_read_support_pass", "high15_final_rank"],
            ascending=[False, True],
            na_position="last",
        ).to_csv(args.outdir / "reference_comparison.tsv", sep="\t", index=False)
    sample_qc.to_csv(args.outdir / "sample_qc.tsv", sep="\t", index=False)
    pd.DataFrame([summary]).to_csv(args.outdir / "run_summary.tsv", sep="\t", index=False)
    (args.outdir / "run_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"Analysis complete: {report}")
    print(f"Main table: {args.outdir / 'utr_results_easy.tsv'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
