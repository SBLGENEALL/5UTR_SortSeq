import sys
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from compare_metrics import analyze_metric_comparison  # noqa: E402


def row(
    variant_id: str,
    total: int,
    unsorted: int,
    high_count: int,
    score: float,
    top15_log2: float,
    high15_probability: float,
    reference: bool = False,
) -> dict[str, object]:
    probabilities = [
        high15_probability / 2,
        high15_probability / 2,
        0.20,
        0.20,
        (0.60 - high15_probability) / 2,
        (0.60 - high15_probability) / 2,
    ]
    result: dict[str, object] = {
        "variant_id": variant_id,
        "total_6bin_count": total,
        "unsorted_count": unsorted,
        "high_bin_raw_count": high_count,
        "expected_bin_score": score,
        "high15_probability": high15_probability,
        "top15_vs_unsorted_log2_enrichment": top15_log2,
        "is_reference_variant": reference,
    }
    for number, probability in enumerate(probabilities, start=1):
        result[f"bin{number}_probability"] = probability
    return result


class CompareMetricsTest(unittest.TestCase):
    def setUp(self):
        self.input = pd.DataFrame(
            [
                row("excluded_200", 200, 100, 50, 6.0, 6.0, 0.60),
                row("orginal", 500, 100, 50, 3.0, 0.0, 0.15, True),
                row("both_high", 500, 100, 50, 5.0, 3.0, 0.45),
                row("score_only", 500, 100, 50, 4.5, -1.0, 0.10),
                row("top15_only", 500, 100, 50, 2.5, 2.0, 0.35),
                row("both_low", 500, 100, 50, 2.0, -2.0, 0.05),
                row("low_unsorted", 500, 10, 50, 5.5, 5.0, 0.55),
            ]
        )

    def test_filter_ranks_overlap_and_reference_quadrants(self):
        all_rows, supported, summary, overlap, summary_json = (
            analyze_metric_comparison(
                self.input,
                min_total_count=201,
                min_unsorted_count=50,
                min_high_bin_count=20,
                top_n=2,
            )
        )
        indexed = all_rows.set_index("variant_id")
        self.assertFalse(indexed.loc["excluded_200", "total_count_filter_pass"])
        self.assertTrue(indexed.loc["orginal", "total_count_filter_pass"])
        self.assertFalse(
            indexed.loc["low_unsorted", "comparison_read_support_pass"]
        )
        self.assertEqual(len(supported), 5)

        supported_indexed = supported.set_index("variant_id")
        self.assertEqual(
            supported_indexed.loc["both_high", "top_list_membership"],
            "consensus",
        )
        self.assertEqual(
            supported_indexed.loc["score_only", "reference_quadrant"],
            "score_only_above",
        )
        self.assertEqual(
            supported_indexed.loc["top15_only", "reference_quadrant"],
            "top15_only_above",
        )
        self.assertEqual(summary_json["actual_top_n"], 2)
        self.assertEqual(summary_json["top_list_consensus_count"], 1)
        self.assertEqual(summary_json["top_list_score_only_count"], 1)
        self.assertEqual(summary_json["top_list_top15_only_count"], 1)
        self.assertIn("metric", summary.columns)
        self.assertIn("requested_top_n", overlap.columns)
        self.assertTrue(summary_json["step8_identity_check_pass"])
        self.assertEqual(summary_json["step6_step8_rank_mismatch_count"], 0)
        self.assertAlmostEqual(
            summary_json["step6_vs_step8_spearman_rho"], 1.0, places=12
        )
        np.testing.assert_allclose(
            supported["step8_high15_relative_enrichment"],
            supported["step6_high15_probability"] / 0.15,
        )
        step6_step8_top2 = overlap[
            (overlap["pair"] == "step6_vs_step8")
            & (overlap["requested_top_n"] == 2)
        ].iloc[0]
        self.assertEqual(step6_step8_top2["overlap_count"], 2)
        self.assertEqual(
            step6_step8_top2["overlap_percent_of_each_list"], 100.0
        )

    def test_missing_metric_fails(self):
        with self.assertRaisesRegex(
            ValueError, "high15_probability"
        ):
            analyze_metric_comparison(
                self.input.drop(columns=["high15_probability"])
            )

    def test_cli_writes_explicit_step678_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            input_path = root / "results.tsv"
            output_dir = root / "comparison"
            self.input.to_csv(input_path, sep="\t", index=False)
            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "compare_metrics.py"),
                    "--input",
                    str(input_path),
                    "--outdir",
                    str(output_dir),
                    "--min-total-count",
                    "201",
                    "--min-unsorted-count",
                    "50",
                    "--min-high-bin-count",
                    "20",
                    "--top-n",
                    "2",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            for filename in [
                "step6_step7_step8_metrics.csv",
                "step6_step7_step8_correlations.csv",
                "step6_step7_step8_topn_overlap.csv",
                "step6_step7_step8_summary.json",
            ]:
                self.assertTrue((output_dir / filename).exists(), filename)


if __name__ == "__main__":
    unittest.main()
