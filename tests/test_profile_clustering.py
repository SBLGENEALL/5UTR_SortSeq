import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "cluster_utr_profiles.py"
SPEC = importlib.util.spec_from_file_location("cluster_utr_profiles", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class ProfileClusteringTests(unittest.TestCase):
    def make_result(self):
        rows = []
        profiles = {
            "high": np.array([0.55, 0.25, 0.10, 0.05, 0.03, 0.02]),
            "middle": np.array([0.05, 0.10, 0.35, 0.30, 0.15, 0.05]),
            "low": np.array([0.01, 0.02, 0.05, 0.12, 0.40, 0.40]),
        }
        weights = np.array([6, 5, 4, 3, 2, 1])
        rng = np.random.default_rng(17)
        for group, profile in profiles.items():
            for number in range(4):
                perturbed = np.clip(profile + rng.normal(0, 0.002, 6), 0, None)
                perturbed = perturbed / perturbed.sum()
                row = {
                    "variant_id": f"{group}_{number}",
                    "top15_read_support_pass": True,
                    "expected_bin_score": float(perturbed @ weights),
                    "high15_probability": float(perturbed[:2].sum()),
                    "is_reference_variant": group == "middle" and number == 0,
                    "candidate_tier": "not_candidate",
                    "recommended_for_cloning": group == "high",
                    "total_6bin_count": 1000,
                    "unsorted_count": 500,
                }
                for bin_number, value in enumerate(perturbed, start=1):
                    row[f"bin{bin_number}_probability"] = value
                rows.append(row)
        # Whole-library profiling must retain low-expression shapes even when
        # they fail the High15-specific bin1+bin2 support rule.
        rows[-1]["top15_read_support_pass"] = False
        excluded = rows[0].copy()
        excluded["variant_id"] = "excluded_low_coverage"
        excluded["top15_read_support_pass"] = False
        excluded["total_6bin_count"] = 100
        rows.append(excluded)
        return pd.DataFrame(rows)

    def test_profiles_are_clustered_and_ordered_high_to_low(self):
        assignments, summary, manifest = MODULE.cluster_profiles(
            self.make_result(), requested_clusters=3
        )
        self.assertEqual(manifest["variants_input"], 13)
        self.assertEqual(manifest["variants_included"], 12)
        self.assertEqual(manifest["actual_clusters"], 3)
        self.assertEqual(manifest["eligibility_rule"], "total_6bin_count >= 200")
        self.assertNotIn("excluded_low_coverage", set(assignments["variant_id"]))
        self.assertIn("low_3", set(assignments["variant_id"]))
        np.testing.assert_allclose(
            assignments[[f"bin{x}_probability" for x in range(1, 7)]].sum(axis=1),
            np.ones(12),
        )
        np.testing.assert_allclose(
            assignments[
                [f"bin{x}_normalized_enrichment_share" for x in range(1, 7)]
            ].sum(axis=1),
            np.ones(12),
        )
        self.assertEqual(sorted(assignments["heatmap_order"]), list(range(1, 13)))
        self.assertGreater(
            summary.iloc[0]["median_equal_bin_high_share"],
            summary.iloc[-1]["median_equal_bin_high_share"],
        )
        high_clusters = set(
            assignments.loc[
                assignments["variant_id"].str.startswith("high_"),
                "profile_cluster",
            ]
        )
        low_clusters = set(
            assignments.loc[
                assignments["variant_id"].str.startswith("low_"),
                "profile_cluster",
            ]
        )
        self.assertEqual(len(high_clusters), 1)
        self.assertEqual(len(low_clusters), 1)
        self.assertNotEqual(high_clusters, low_clusters)
        reference = assignments[assignments["is_reference_variant"]]
        self.assertEqual(reference["variant_id"].tolist(), ["middle_0"])

    def test_all_utrs_are_exported_and_top_candidates_include_reference(self):
        assignments, summary, all_profiles, top_profiles, manifest = (
            MODULE.build_profile_outputs(
                self.make_result(), requested_clusters=3, top_n=2
            )
        )
        self.assertEqual(len(assignments), 12)
        self.assertEqual(len(summary), 3)
        self.assertEqual(len(all_profiles), 13)
        self.assertEqual(manifest["variants_exported_all"], 13)
        self.assertEqual(int(top_profiles["top_candidate_selected"].sum()), 2)
        reference = top_profiles[top_profiles["is_reference_variant"]]
        self.assertEqual(reference["variant_id"].tolist(), ["middle_0"])
        self.assertTrue(reference["reference_added_for_plot"].iloc[0])
        q_columns = [
            f"bin{x}_normalized_enrichment_share" for x in range(1, 7)
        ]
        np.testing.assert_allclose(
            all_profiles[q_columns].sum(axis=1), np.ones(len(all_profiles))
        )

    def test_invalid_cluster_count_fails(self):
        with self.assertRaisesRegex(ValueError, "at least 2"):
            MODULE.cluster_profiles(self.make_result(), requested_clusters=1)

    def test_invalid_total_count_cutoff_fails(self):
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            MODULE.cluster_profiles(
                self.make_result(), requested_clusters=3, min_total_count=-1
            )


if __name__ == "__main__":
    unittest.main()
