import sys
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from qc_bimodality import analyze_bimodality  # noqa: E402


class BimodalityQCTest(unittest.TestCase):
    def setUp(self):
        self.input = pd.DataFrame(
            [
                {
                    "variant_id": "UTR_A",
                    "total_6bin_count": 1500,
                    "unsorted_count": 500,
                    "high_bin_raw_count": 300,
                    "bin1_count": 150,
                    "bin2_count": 150,
                    "bin5_count": 300,
                    "bin6_count": 450,
                    "bin1_probability": 0.20,
                    "bin2_probability": 0.20,
                    "bin3_probability": 0.04,
                    "bin4_probability": 0.06,
                    "bin5_probability": 0.20,
                    "bin6_probability": 0.30,
                },
                {
                    "variant_id": "middle",
                    "total_6bin_count": 1500,
                    "unsorted_count": 500,
                    "high_bin_raw_count": 300,
                    "bin1_count": 150,
                    "bin2_count": 150,
                    "bin5_count": 0,
                    "bin6_count": 0,
                    "bin1_probability": 0.05,
                    "bin2_probability": 0.10,
                    "bin3_probability": 0.60,
                    "bin4_probability": 0.25,
                    "bin5_probability": 0.00,
                    "bin6_probability": 0.00,
                },
                {
                    "variant_id": "neutral",
                    "total_6bin_count": 1500,
                    "unsorted_count": 500,
                    "high_bin_raw_count": 225,
                    "bin1_count": 75,
                    "bin2_count": 150,
                    "bin5_count": 450,
                    "bin6_count": 300,
                    "bin1_probability": 0.05,
                    "bin2_probability": 0.10,
                    "bin3_probability": 0.15,
                    "bin4_probability": 0.20,
                    "bin5_probability": 0.30,
                    "bin6_probability": 0.20,
                },
                {
                    "variant_id": "sparse_polarized",
                    "total_6bin_count": 100,
                    "unsorted_count": 20,
                    "high_bin_raw_count": 10,
                    "bin1_count": 5,
                    "bin2_count": 5,
                    "bin5_count": 20,
                    "bin6_count": 30,
                    "bin1_probability": 0.20,
                    "bin2_probability": 0.20,
                    "bin3_probability": 0.04,
                    "bin4_probability": 0.06,
                    "bin5_probability": 0.20,
                    "bin6_probability": 0.30,
                },
            ]
        )

    def test_shape_and_read_support_are_separate(self):
        result, by_count, sensitivity, summary = analyze_bimodality(
            self.input,
            min_total_count=200,
            min_each_tail_count=20,
            min_high_tail_probability=0.20,
            min_low_tail_probability=0.20,
            max_middle_probability=0.30,
            max_valley_ratio=0.75,
        )
        indexed = result.set_index("variant_id")

        self.assertTrue(indexed.loc["UTR_A", "clear_bimodal_flag"])
        self.assertTrue(
            indexed.loc["UTR_A", "utra_like_strong_polarization_flag"]
        )
        self.assertAlmostEqual(
            indexed.loc["UTR_A", "extreme_polarization_index_0to1"], 0.8
        )
        self.assertFalse(indexed.loc["middle", "clear_bimodal_flag"])
        self.assertFalse(indexed.loc["middle", "both_tail_raw_support_pass"])
        self.assertFalse(indexed.loc["neutral", "clear_bimodal_flag"])

        self.assertTrue(
            indexed.loc["sparse_polarized", "extreme_polarized_shape_flag"]
        )
        self.assertFalse(
            indexed.loc["sparse_polarized", "bimodality_read_support_pass"]
        )
        self.assertFalse(
            indexed.loc["sparse_polarized", "clear_bimodal_flag"]
        )

        self.assertEqual(summary["variants_total"], 4)
        self.assertEqual(summary["variants_read_supported"], 3)
        self.assertEqual(summary["clear_bimodal_count"], 1)
        self.assertEqual(summary["utra_like_strong_polarization_count"], 1)
        self.assertIn("clear_bimodal_percent", sensitivity.columns)
        self.assertEqual(
            int(by_count.loc[by_count["read_count_group"] == "<200", "all_variants"].iloc[0]),
            1,
        )

    def test_missing_probability_column_fails(self):
        with self.assertRaisesRegex(ValueError, "bin6_probability"):
            analyze_bimodality(
                self.input.drop(columns=["bin6_probability"]),
                200,
                20,
                0.20,
                0.20,
                0.30,
                0.75,
            )


if __name__ == "__main__":
    unittest.main()
