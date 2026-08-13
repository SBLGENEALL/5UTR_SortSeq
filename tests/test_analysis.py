import importlib.util
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
        self.assertAlmostEqual(neutral["high15_enrichment"], 1.0, places=6)

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

    def test_collected_cells_override_nominal_population_fractions(self):
        _, samples, _ = self.make_inputs()
        samples["count_file"] = "unused.tsv"
        samples["cells_collected"] = [np.nan, 10, 10, 10, 10, 10, 50]
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "samples.tsv"
            samples.to_csv(path, sep="\t", index=False)
            loaded = MODULE.load_samples(path)
        bins = loaded[loaded["sample_type"].eq("bin")].sort_values("bin_number")
        np.testing.assert_allclose(bins["population_fraction"], [0.10] * 5 + [0.50])
        self.assertEqual(set(bins["population_fraction_source"]), {"cells_collected"})


if __name__ == "__main__":
    unittest.main()
