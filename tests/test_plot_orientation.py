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
            "29_top_high15_individual_relative_enrichment_profiles.png", script
        )


if __name__ == "__main__":
    unittest.main()
