#!/usr/bin/env python3
"""Quantify U-shaped/high+low-tail Sort-seq distributions and count dependence."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="QC high+low-tail polarization and its dependence on read support."
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--min-total-count", type=int, default=200)
    parser.add_argument("--min-each-tail-count", type=int, default=20)
    parser.add_argument("--min-high-tail-probability", type=float, default=0.20)
    parser.add_argument("--min-low-tail-probability", type=float, default=0.20)
    parser.add_argument("--max-middle-probability", type=float, default=0.30)
    parser.add_argument("--max-valley-ratio", type=float, default=0.75)
    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    separator = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=separator)


def write_csv(table: pd.DataFrame, path: Path) -> None:
    table.to_csv(path, index=False, encoding="utf-8-sig")


def count_group(values: pd.Series) -> pd.Categorical:
    return pd.cut(
        values,
        bins=[-np.inf, 199, 499, 999, 4999, np.inf],
        labels=["<200", "200-499", "500-999", "1000-4999", ">=5000"],
        ordered=True,
    )


def analyze_bimodality(
    result: pd.DataFrame,
    min_total_count: int,
    min_each_tail_count: int,
    min_high_tail_probability: float,
    min_low_tail_probability: float,
    max_middle_probability: float,
    max_valley_ratio: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    probability_columns = [f"bin{x}_probability" for x in range(1, 7)]
    required = {
        "variant_id",
        "total_6bin_count",
        "unsorted_count",
        "bin1_count",
        "bin2_count",
        "bin5_count",
        "bin6_count",
        *probability_columns,
    }
    missing = sorted(required.difference(result.columns))
    if missing:
        raise ValueError(f"input is missing columns: {missing}")
    if min_total_count < 0 or min_each_tail_count < 0:
        raise ValueError("read-count cutoffs cannot be negative")
    for value, label in [
        (min_high_tail_probability, "min_high_tail_probability"),
        (min_low_tail_probability, "min_low_tail_probability"),
        (max_middle_probability, "max_middle_probability"),
    ]:
        if not 0 <= value <= 1:
            raise ValueError(f"{label} must be between 0 and 1")
    if max_valley_ratio < 0:
        raise ValueError("max_valley_ratio cannot be negative")

    output = result.copy()
    for column in [
        "total_6bin_count",
        "unsorted_count",
        "bin1_count",
        "bin2_count",
        "bin5_count",
        "bin6_count",
        *probability_columns,
    ]:
        output[column] = pd.to_numeric(output[column], errors="coerce")

    p = output[probability_columns]
    output["high_tail_probability"] = p.iloc[:, :2].sum(axis=1, min_count=2)
    output["middle_probability"] = p.iloc[:, 2:4].sum(axis=1, min_count=2)
    output["low_tail_probability"] = p.iloc[:, 4:].sum(axis=1, min_count=2)
    output["extreme_tail_mass"] = (
        output["high_tail_probability"] + output["low_tail_probability"]
    )
    tail_sum = output["extreme_tail_mass"].replace(0, np.nan)
    output["tail_balance_0to1"] = (
        2
        * np.minimum(
            output["high_tail_probability"], output["low_tail_probability"]
        )
        / tail_sum
    )
    output["extreme_polarization_index_0to1"] = (
        output["extreme_tail_mass"] * output["tail_balance_0to1"]
    )

    output["high_tail_peak"] = p.iloc[:, :2].max(axis=1)
    output["middle_peak"] = p.iloc[:, 2:4].max(axis=1)
    output["low_tail_peak"] = p.iloc[:, 4:].max(axis=1)
    smaller_tail_peak = np.minimum(
        output["high_tail_peak"], output["low_tail_peak"]
    ).replace(0, np.nan)
    output["middle_valley_ratio"] = output["middle_peak"] / smaller_tail_peak
    output["high_tail_raw_count"] = output["bin1_count"] + output["bin2_count"]
    output["low_tail_raw_count"] = output["bin5_count"] + output["bin6_count"]

    finite_distribution = p.notna().all(axis=1)
    output["bimodality_read_support_pass"] = (
        finite_distribution
        & (output["total_6bin_count"] >= min_total_count)
    )
    output["both_tail_raw_support_pass"] = (
        finite_distribution
        & (output["high_tail_raw_count"] >= min_each_tail_count)
        & (output["low_tail_raw_count"] >= min_each_tail_count)
    )
    output["extreme_polarized_shape_flag"] = (
        finite_distribution
        & (output["high_tail_probability"] >= min_high_tail_probability)
        & (output["low_tail_probability"] >= min_low_tail_probability)
        & (output["middle_probability"] <= max_middle_probability)
    )
    output["clear_bimodal_flag"] = (
        output["bimodality_read_support_pass"]
        & output["both_tail_raw_support_pass"]
        & output["extreme_polarized_shape_flag"]
        & (output["middle_valley_ratio"] <= max_valley_ratio)
    )
    output["utra_like_strong_polarization_flag"] = (
        output["bimodality_read_support_pass"]
        & output["both_tail_raw_support_pass"]
        & (output["high_tail_probability"] >= 0.20)
        & (output["low_tail_probability"] >= 0.50)
        & (output["middle_probability"] <= 0.20)
        & (output["middle_valley_ratio"] <= max_valley_ratio)
    )
    output["read_count_group"] = count_group(output["total_6bin_count"])

    group_rows: list[dict[str, object]] = []
    for group in output["read_count_group"].cat.categories:
        subset = output[output["read_count_group"].eq(group)]
        supported = subset[subset["bimodality_read_support_pass"]]
        group_rows.append(
            {
                "read_count_group": str(group),
                "all_variants": int(len(subset)),
                "read_supported_variants": int(len(supported)),
                "clear_bimodal_count": int(supported["clear_bimodal_flag"].sum()),
                "clear_bimodal_percent_of_supported": (
                    100 * float(supported["clear_bimodal_flag"].mean())
                    if len(supported)
                    else np.nan
                ),
                "utra_like_count": int(
                    supported["utra_like_strong_polarization_flag"].sum()
                ),
                "median_polarization_index": (
                    float(supported["extreme_polarization_index_0to1"].median())
                    if len(supported)
                    else np.nan
                ),
            }
        )
    by_count = pd.DataFrame(group_rows)

    sensitivity_rows: list[dict[str, object]] = []
    for threshold in [100, 200, 500, 1000, 2000, 5000]:
        subset = output[
            finite_distribution
            & (output["total_6bin_count"] >= threshold)
        ]
        shape = (
            subset["both_tail_raw_support_pass"]
            & subset["extreme_polarized_shape_flag"]
            & (subset["middle_valley_ratio"] <= max_valley_ratio)
        )
        sensitivity_rows.append(
            {
                "minimum_total_6bin_count": threshold,
                "variants_passing": int(len(subset)),
                "clear_bimodal_count": int(shape.sum()),
                "clear_bimodal_percent": (
                    100 * float(shape.mean()) if len(subset) else np.nan
                ),
            }
        )
    sensitivity = pd.DataFrame(sensitivity_rows)

    supported = output[output["bimodality_read_support_pass"]]
    if (
        len(supported) >= 3
        and supported["total_6bin_count"].nunique() >= 2
        and supported["extreme_polarization_index_0to1"].nunique() >= 2
    ):
        correlation = spearmanr(
            np.log10(supported["total_6bin_count"] + 1),
            supported["extreme_polarization_index_0to1"],
            nan_policy="omit",
        )
        rho = float(correlation.statistic)
        pvalue = float(correlation.pvalue)
    else:
        rho = float("nan")
        pvalue = float("nan")

    summary = {
        "variants_total": int(len(output)),
        "variants_read_supported": int(output["bimodality_read_support_pass"].sum()),
        "extreme_polarized_shape_count_all_coverage": int(
            output["extreme_polarized_shape_flag"].sum()
        ),
        "clear_bimodal_count": int(output["clear_bimodal_flag"].sum()),
        "clear_bimodal_percent_of_supported": (
            100 * float(supported["clear_bimodal_flag"].mean())
            if len(supported)
            else float("nan")
        ),
        "utra_like_strong_polarization_count": int(
            output["utra_like_strong_polarization_flag"].sum()
        ),
        "spearman_log10_total_count_vs_polarization_index": rho,
        "spearman_pvalue": pvalue,
        "min_total_count": min_total_count,
        "min_each_tail_count": min_each_tail_count,
        "min_high_tail_probability": min_high_tail_probability,
        "min_low_tail_probability": min_low_tail_probability,
        "max_middle_probability": max_middle_probability,
        "max_valley_ratio": max_valley_ratio,
        "interpretation": (
            "Clear bimodality is an exploratory shape/QC flag, not proof of a "
            "UTR-intrinsic biological state."
        ),
    }
    output = output.sort_values(
        ["clear_bimodal_flag", "extreme_polarization_index_0to1"],
        ascending=[False, False],
        na_position="last",
    )
    return output, by_count, sensitivity, summary


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result = read_table(args.input)
    output, by_count, sensitivity, summary = analyze_bimodality(
        result,
        args.min_total_count,
        args.min_each_tail_count,
        args.min_high_tail_probability,
        args.min_low_tail_probability,
        args.max_middle_probability,
        args.max_valley_ratio,
    )
    write_csv(output, args.outdir / "bimodality_all_variants.csv")
    write_csv(
        output[output["clear_bimodal_flag"]],
        args.outdir / "clear_bimodal_candidates.csv",
    )
    write_csv(
        output[output["utra_like_strong_polarization_flag"]],
        args.outdir / "utra_like_strong_polarization.csv",
    )
    write_csv(by_count, args.outdir / "bimodality_by_read_count.csv")
    write_csv(sensitivity, args.outdir / "bimodality_count_sensitivity.csv")
    write_csv(pd.DataFrame([summary]), args.outdir / "bimodality_summary.csv")
    (args.outdir / "bimodality_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"Bimodality QC written to: {args.outdir}")
    print(
        "Clear bimodal: "
        f"{summary['clear_bimodal_count']:,} / "
        f"{summary['variants_read_supported']:,} supported "
        f"({summary['clear_bimodal_percent_of_supported']:.2f}%)"
    )
    print(
        "UTR-A-like strong polarization: "
        f"{summary['utra_like_strong_polarization_count']:,}"
    )
    print(
        "Read-count correlation (Spearman rho): "
        f"{summary['spearman_log10_total_count_vs_polarization_index']:.4f}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
