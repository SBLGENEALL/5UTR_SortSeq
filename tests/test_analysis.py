import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "analyze_sortseq.py"
SPEC = importlib.util.spec_from_file_location("analyze_sortseq", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class EasyPipelineTests(unittest.TestCase):
    def make_inputs(self):
        variants = pd.DataFrame({"variant_id": ["high", "neutral", "low"]})
        samples = pd.DataFrame(
            {
                "sample_id": ["unsorted", "bin1", "bin2", "bin3", "bin4", "bin5", "bin6"],
                "sample_type": ["unsorted"] + ["bin"] * 6,
                "bin_number": [np.nan, 1, 2, 3, 4, 5, 6],
                "population_fraction": [np.nan, 0.05, 0.10, 0.15, 0.20, 0.30, 0.20],
                "representative_log10_mfi": [np.nan, 6, 5, 4, 3, 2, 1],
            }
        )
        # Different total depths in every bin.  Within-bin composition is set so
        # high shifts toward bin1, neutral stays at pool frequency, and low shifts
        # toward bin6.
        frequencies = {
            "bin1": [0.70, 0.20, 0.10],
            "bin2": [0.60, 0.20, 0.20],
            "bin3": [0.45, 0.20, 0.35],
            "bin4": [0.30, 0.20, 0.50],
            "bin5": [0.20, 0.20, 0.60],
            "bin6": [0.10, 0.20, 0.70],
        }
        depths = {
            "bin1": 1000,
            "bin2": 2000,
            "bin3": 3000,
            "bin4": 4000,
            "bin5": 5000,
            "bin6": 6000,
        }
        counts = pd.DataFrame(index=pd.Index(["high", "neutral", "low"], name="variant_id"))
        # Neutral is 20% in both whole unsorted and every bin, so its estimated
        # target-gate entry probability should equal the global 90% gate rate.
        counts["unsorted"] = [400, 200, 400]
        for sample_id, frequency in frequencies.items():
            counts[sample_id] = np.asarray(frequency) * depths[sample_id]
        return variants, samples, counts

    def test_high_neutral_low_order_is_recovered(self):
        variants, samples, counts = self.make_inputs()
        result, _, summary = MODULE.analyze(
            variants,
            samples,
            counts,
            "variant_id",
            0.5,
            0,
            0,
            0.90,
        )
        result = result.set_index("variant_id")
        self.assertGreater(
            result.loc["high", "expected_bin_score"],
            result.loc["neutral", "expected_bin_score"],
        )
        self.assertGreater(
            result.loc["neutral", "expected_bin_score"],
            result.loc["low", "expected_bin_score"],
        )
        self.assertAlmostEqual(result.loc["neutral", "expected_bin_score"], 2.8, places=6)
        self.assertAlmostEqual(summary["neutral_baseline_score"], 2.8, places=6)

    def test_neutral_distribution_matches_population_fractions(self):
        variants, samples, counts = self.make_inputs()
        result, _, _ = MODULE.analyze(
            variants,
            samples,
            counts,
            "variant_id",
            0.5,
            0,
            0,
            0.90,
        )
        neutral = result.set_index("variant_id").loc["neutral"]
        expected = [0.05, 0.10, 0.15, 0.20, 0.30, 0.20]
        observed = [neutral[f"bin{x}_probability"] for x in range(1, 7)]
        np.testing.assert_allclose(observed, expected)
        normalized = [
            neutral[f"bin{x}_normalized_enrichment_share"]
            for x in range(1, 7)
        ]
        np.testing.assert_allclose(normalized, np.repeat(1 / 6, 6))
        self.assertAlmostEqual(
            neutral["normalized_enrichment_share_sum"], 1.0, places=12
        )
        self.assertAlmostEqual(neutral["equal_bin_high_share"], 1 / 3, places=12)
        self.assertAlmostEqual(neutral["high15_enrichment"], 1.0, places=6)
        self.assertAlmostEqual(
            neutral["top15_vs_unsorted_enrichment"], 1.0, delta=0.003
        )
        self.assertGreater(
            result.set_index("variant_id").loc[
                "high", "top15_vs_unsorted_enrichment"
            ],
            neutral["top15_vs_unsorted_enrichment"],
        )

    def test_whole_unsorted_gate_probability_is_secondary_estimate(self):
        variants, samples, counts = self.make_inputs()
        result, _, _ = MODULE.analyze(
            variants,
            samples,
            counts,
            "variant_id",
            0.5,
            0,
            0,
            0.90,
        )
        neutral = result.set_index("variant_id").loc["neutral"]
        self.assertAlmostEqual(neutral["gate_entry_probability_raw"], 0.90, delta=0.001)

    def test_strict_filter_separates_sparse_from_supported_candidates(self):
        variants, samples, counts = self.make_inputs()
        result, _, summary = MODULE.analyze(
            variants,
            samples,
            counts,
            "variant_id",
            0.5,
            0,
            0,
            0.90,
            strict_min_unsorted_count=300,
            strict_min_total_bin_count=5000,
            strict_relative_median_fraction=0.0,
            strict_min_high_bin_count=200,
            strict_min_detected_bins=3,
        )
        indexed = result.set_index("variant_id")
        self.assertTrue(indexed.loc["high", "strict_coverage_pass"])
        self.assertFalse(indexed.loc["neutral", "strict_coverage_pass"])
        self.assertTrue(indexed.loc["high", "high_confidence_candidate_flag"])
        self.assertEqual(summary["strict_unsorted_cutoff"], 300)
        self.assertEqual(summary["strict_total_6bin_cutoff"], 5000)
        self.assertGreaterEqual(summary["strict_coverage_passing_count"], 1)

    def test_orginal_alias_adds_reference_relative_metrics(self):
        variants, samples, counts = self.make_inputs()
        variants.loc[variants["variant_id"].eq("neutral"), "variant_id"] = "orginal"
        counts = counts.rename(index={"neutral": "orginal"})
        reference_id = MODULE.resolve_reference_variant_id(variants, "variant_id", "auto")
        self.assertEqual(reference_id, "orginal")
        result, essential, summary = MODULE.analyze(
            variants,
            samples,
            counts,
            "variant_id",
            0.5,
            0,
            0,
            0.90,
            reference_id,
        )
        indexed = result.set_index("variant_id")
        self.assertTrue(indexed.loc["orginal", "is_reference_variant"])
        self.assertAlmostEqual(indexed.loc["orginal", "delta_score_vs_reference"], 0.0)
        self.assertGreater(indexed.loc["high", "delta_score_vs_reference"], 0)
        self.assertLess(indexed.loc["low", "delta_score_vs_reference"], 0)
        self.assertGreater(indexed.loc["high", "high15_fold_vs_reference"], 1)
        self.assertGreater(
            indexed.loc["high", "delta_top15_log2_enrichment_vs_reference"], 0
        )
        self.assertTrue(indexed.loc["high", "top15_candidate_flag"])
        self.assertAlmostEqual(
            indexed.loc["orginal", "delta_top15_log2_enrichment_vs_reference"],
            0.0,
            places=10,
        )
        self.assertEqual(summary["reference_variant_id"], "orginal")
        self.assertEqual(summary["variants_with_score_above_reference"], 1)
        self.assertIn("delta_score_vs_reference", essential.columns)

    def test_reference_alias_can_be_found_in_collapsed_original_ids(self):
        variants = pd.DataFrame(
            {
                "variant_id": ["sequence_key_1", "sequence_key_2"],
                "original_variant_ids": ["control|orginal", "other"],
            }
        )
        reference_id = MODULE.resolve_reference_variant_id(variants, "variant_id", "auto")
        self.assertEqual(reference_id, "sequence_key_1")
        self.assertEqual(
            MODULE.reference_display_label(variants, "variant_id", reference_id),
            "orginal",
        )

    def test_nominal_population_fractions_take_priority_over_collected_events(self):
        _, samples, _ = self.make_inputs()
        samples["count_file"] = "unused.tsv"
        samples["cells_collected"] = [np.nan, 10, 10, 10, 10, 10, 50]
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "samples.tsv"
            samples.to_csv(path, sep="\t", index=False)
            loaded = MODULE.load_samples(path)
        bins = loaded[loaded["sample_type"].eq("bin")].sort_values("bin_number")
        np.testing.assert_allclose(
            bins["population_fraction"], [0.05, 0.10, 0.15, 0.20, 0.30, 0.20]
        )
        np.testing.assert_allclose(
            bins["cells_collected_fraction"], [0.10] * 5 + [0.50]
        )
        self.assertEqual(
            set(bins["population_fraction_source"]), {"population_fraction"}
        )

    def test_high15_primary_is_invariant_to_unsorted_abundance(self):
        variants = pd.DataFrame({"variant_id": ["A", "B", "original"]})
        _, samples, _ = self.make_inputs()
        counts = pd.DataFrame(
            {
                "unsorted": [10, 1000, 100],
                **{f"bin{x}": [100, 100, 100] for x in range(1, 7)},
            },
            index=pd.Index(["A", "B", "original"], name="variant_id"),
        )
        result, _, _ = MODULE.analyze(
            variants,
            samples,
            counts,
            "variant_id",
            0.5,
            0,
            0,
            0.90,
            reference_variant_id="original",
            top_hit_min_unsorted_count=0,
            top_hit_min_total_bin_count=0,
            top_hit_min_high_bin_count=0,
        )
        indexed = result.set_index("variant_id")
        self.assertAlmostEqual(
            indexed.loc["A", "high15_probability"],
            indexed.loc["B", "high15_probability"],
        )
        self.assertGreater(
            indexed.loc["A", "top15_vs_unsorted_enrichment"],
            indexed.loc["B", "top15_vs_unsorted_enrichment"],
        )

    def test_bin1_concentrated_hit_does_not_require_bin2_enrichment(self):
        variants = pd.DataFrame({"variant_id": ["high", "original", "other"]})
        _, samples, _ = self.make_inputs()
        frequencies = {
            "bin1": [0.90, 0.05, 0.05],
            "bin2": [0.01, 0.495, 0.495],
            "bin3": [0.10, 0.45, 0.45],
            "bin4": [0.10, 0.45, 0.45],
            "bin5": [0.10, 0.45, 0.45],
            "bin6": [0.10, 0.45, 0.45],
        }
        counts = pd.DataFrame(
            {"unsorted": [1000, 1000, 1000]},
            index=pd.Index(["high", "original", "other"], name="variant_id"),
        )
        for sample_id, frequency in frequencies.items():
            counts[sample_id] = np.asarray(frequency) * 10000
        result, _, _ = MODULE.analyze(
            variants,
            samples,
            counts,
            "variant_id",
            0.5,
            0,
            0,
            0.90,
            reference_variant_id="original",
            top_hit_min_unsorted_count=0,
            top_hit_min_total_bin_count=0,
            top_hit_min_high_bin_count=0,
        )
        high = result.set_index("variant_id").loc["high"]
        self.assertFalse(high["top15_both_bins_enriched"])
        self.assertTrue(high["top15_above_comparator"])
        self.assertTrue(high["top15_candidate_flag"])

    def test_technical_bootstrap_is_deterministic_and_reference_aware(self):
        raw_bins = pd.DataFrame(
            {
                "bin1": [900, 50, 50],
                "bin2": [100, 450, 450],
                "bin3": [100, 450, 450],
                "bin4": [100, 450, 450],
                "bin5": [100, 450, 450],
                "bin6": [100, 450, 450],
            },
            index=pd.Index(["high", "original", "other"], name="variant_id"),
        )
        fractions = pd.Series(
            [0.05, 0.10, 0.15, 0.20, 0.30, 0.20], index=raw_bins.columns
        )
        support = pd.Series(True, index=raw_bins.index)
        kwargs = dict(
            raw_bins=raw_bins,
            fractions=fractions,
            reference_variant_id="original",
            read_support=support,
            replicates=200,
            seed=1234,
            lower_quantile=0.10,
            top_n=1,
        )
        first = MODULE.technical_high15_bootstrap(**kwargs)
        second = MODULE.technical_high15_bootstrap(**kwargs)
        pd.testing.assert_frame_equal(first, second)
        self.assertGreater(
            first.loc["high", "high15_bootstrap_probability_above_comparator"],
            0.90,
        )
        self.assertGreater(
            first.loc["high", "high15_bootstrap_top_n_frequency"], 0.90
        )

    def test_cli_writes_primary_ranking_and_cloning_tables(self):
        variants, samples, counts = self.make_inputs()
        variants.loc[variants["variant_id"].eq("neutral"), "variant_id"] = "original"
        counts = counts.rename(index={"neutral": "original"})
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            variants_path = root / "variants.tsv"
            samples_path = root / "samples.tsv"
            outdir = root / "result"
            variants.to_csv(variants_path, sep="\t", index=False)
            count_paths = []
            for sample_id in samples["sample_id"]:
                count_path = root / f"{sample_id}.tsv"
                pd.DataFrame(
                    {
                        "variant_id": counts.index,
                        "count": counts[sample_id].to_numpy(),
                    }
                ).to_csv(count_path, sep="\t", index=False)
                count_paths.append(str(count_path))
            samples = samples.copy()
            samples["count_file"] = count_paths
            samples.to_csv(samples_path, sep="\t", index=False)
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--variants",
                    str(variants_path),
                    "--samples",
                    str(samples_path),
                    "--outdir",
                    str(outdir),
                    "--reference-variant-id",
                    "original",
                    "--min-unsorted-count",
                    "0",
                    "--min-total-bin-count",
                    "0",
                    "--top-hit-min-unsorted-count",
                    "0",
                    "--top-hit-min-total-bin-count",
                    "0",
                    "--top-hit-min-high-bin-count",
                    "0",
                    "--bootstrap-replicates",
                    "50",
                    "--bootstrap-top-n",
                    "1",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            ranking = pd.read_csv(outdir / "high15_primary_ranking.csv")
            self.assertEqual(ranking.iloc[0]["variant_id"], "high")
            self.assertIn("high15_robust_rank_score", ranking.columns)
            self.assertIn("bin1_normalized_enrichment_share", ranking.columns)
            self.assertIn("candidate_tier", ranking.columns)
            self.assertTrue((outdir / "top_candidates_for_cloning.csv").exists())
            self.assertTrue((outdir / "top1_candidates_for_cloning.csv").exists())
            legacy = pd.read_csv(outdir / "top15_enrichment_ranking.csv")
            pd.testing.assert_frame_equal(ranking, legacy)


if __name__ == "__main__":
    unittest.main()
