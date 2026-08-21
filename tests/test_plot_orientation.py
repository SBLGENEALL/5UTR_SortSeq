import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class PlotOrientationTests(unittest.TestCase):
    def test_metric_profiles_show_bin1_on_left(self):
        script = (ROOT / "scripts" / "plot_metric_comparison.R").read_text(
            encoding="utf-8"
        )
        self.assertNotIn('levels = paste0("bin", 6:1)', script)
        self.assertIn('levels = paste0("bin", 1:6)', script)
        self.assertIn(
            "29_top200_normalized_enrichment_profile_heatmap.png", script
        )

    def test_normalized_profile_plots_include_top200_and_reference(self):
        script = (ROOT / "scripts" / "plot_sortseq.R").read_text(
            encoding="utf-8"
        )
        for filename in [
            "29_all_utr_normalized_enrichment_profile_heatmap.png",
            "30_normalized_enrichment_profile_cluster_means.png",
            "31_normalized_enrichment_profile_cluster_variability.png",
            "32_top200_normalized_enrichment_profile_heatmap.png",
            "top200_normalized_enrichment_profile_pages.pdf",
        ]:
            self.assertIn(filename, script)
        self.assertIn("dashed red = original/orginal", script)


if __name__ == "__main__":
    unittest.main()
