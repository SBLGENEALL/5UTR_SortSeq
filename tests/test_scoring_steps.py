import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "export_scoring_steps.py"


class ScoringStepExportTests(unittest.TestCase):
    def test_five_steps_recover_neutral_population_distribution(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            variants = ["high", "orginal", "low"]
            frequencies = {
                "B1": [0.70, 0.20, 0.10],
                "B2": [0.60, 0.20, 0.20],
                "B3": [0.45, 0.20, 0.35],
                "B4": [0.30, 0.20, 0.50],
                "B5": [0.20, 0.20, 0.60],
                "B6": [0.10, 0.20, 0.70],
            }
            depths = {"B1": 1000, "B2": 2000, "B3": 3000, "B4": 4000, "B5": 5000, "B6": 6000}
            matrix = pd.DataFrame(
                {
                    "sequence_key": variants,
                    "variant_ids": variants,
                    "U": [400, 200, 400],
                }
            )
            for column, frequency in frequencies.items():
                matrix[column] = (np.asarray(frequency) * depths[column]).astype(int)
            matrix_path = root / "matrix.csv"
            matrix.to_csv(matrix_path, index=False)

            mapping = pd.DataFrame(
                {
                    "sample_id": ["bin1", "bin2", "bin3", "bin4", "bin5", "bin6", "unsorted"],
                    "ngs_column": ["B1", "B2", "B3", "B4", "B5", "B6", "U"],
                    "sample_type": ["bin"] * 6 + ["unsorted"],
                    "bin_number": [1, 2, 3, 4, 5, 6, ""],
                    "population_fraction": [0.05, 0.10, 0.15, 0.20, 0.30, 0.20, ""],
                }
            )
            map_path = root / "sample_map.csv"
            mapping.to_csv(map_path, index=False)
            output = root / "steps"

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--matrix",
                    str(matrix_path),
                    "--sample-map",
                    str(map_path),
                    "--outdir",
                    str(output),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            for number in range(1, 6):
                self.assertTrue(list(output.glob(f"{number:02d}_*.csv")))
            step2 = pd.read_csv(output / "02_depth_normalized_frequency.csv")
            frequency_columns = [f"bin{x}_depth_normalized_frequency" for x in range(1, 7)]
            np.testing.assert_allclose(step2[frequency_columns].sum(axis=0), np.ones(6))

            step4 = pd.read_csv(output / "04_within_utr_bin_probability.csv").set_index("variant_id")
            probability_columns = [f"bin{x}_probability" for x in range(1, 7)]
            np.testing.assert_allclose(step4[probability_columns].sum(axis=1), np.ones(3))
            np.testing.assert_allclose(
                step4.loc["orginal", probability_columns].to_numpy(dtype=float),
                [0.05, 0.10, 0.15, 0.20, 0.30, 0.20],
            )

            step5 = pd.read_csv(output / "05_score_contributions_and_final_score.csv").set_index("variant_id")
            self.assertAlmostEqual(step5.loc["orginal", "expected_bin_score"], 2.8, places=6)
            self.assertAlmostEqual(step5.loc["orginal", "delta_score_vs_reference"], 0.0, places=6)
            self.assertGreater(
                step5.loc["high", "expected_bin_score"],
                step5.loc["orginal", "expected_bin_score"],
            )
            self.assertGreater(
                step5.loc["high", "top15_vs_unsorted_enrichment"],
                step5.loc["orginal", "top15_vs_unsorted_enrichment"],
            )
            self.assertAlmostEqual(
                step5.loc["orginal", "delta_top15_log2_enrichment_vs_reference"],
                0.0,
                places=10,
            )
            ranking = pd.read_csv(
                output / "08_top15_unsorted_enrichment_ranking.csv"
            )
            self.assertEqual(ranking.iloc[0]["variant_id"], "high")


if __name__ == "__main__":
    unittest.main()
