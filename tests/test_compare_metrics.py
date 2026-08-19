import sys
import unittest
from pathlib import Path

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

    def test_missing_metric_fails(self):
        with self.assertRaisesRegex(
            ValueError, "high15_probability"
        ):
            analyze_metric_comparison(
                self.input.drop(columns=["high15_probability"])
            )


if __name__ == "__main__":
    unittest.main()
