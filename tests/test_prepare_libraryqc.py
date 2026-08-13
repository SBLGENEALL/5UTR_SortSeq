import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "prepare_libraryqc_sortseq.py"


class PrepareLibraryQCTests(unittest.TestCase):
    def test_matrix_conversion(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            matrix = pd.DataFrame(
                {
                    "sequence_key": ["UTR_A", "UTR_B"],
                    "variant_ids": ["UTR_A", "UTR_B"],
                    "target_sequence": ["AAAA", "CCCC"],
                    "length": [4, 4],
                    "gc_percent": [0, 100],
                    "B1": [10, 20],
                    "B2": [11, 21],
                    "B3": [12, 22],
                    "B4": [13, 23],
                    "B5": [14, 24],
                    "B6": [15, 25],
                    "U": [16, 26],
                    "Undetermined_residual": [99, 99],
                }
            )
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
            map_path = root / "map.csv"
            mapping.to_csv(map_path, index=False)
            outdir = root / "out"
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--matrix",
                    str(matrix_path),
                    "--sample-map",
                    str(map_path),
                    "--outdir",
                    str(outdir),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            bin1 = pd.read_csv(outdir / "counts" / "bin1.tsv", sep="\t")
            self.assertEqual(bin1["count"].tolist(), [10, 20])
            samples = pd.read_csv(outdir / "samples.tsv", sep="\t")
            self.assertEqual(len(samples), 7)
            self.assertNotIn("Undetermined_residual", samples["sample_id"].tolist())


if __name__ == "__main__":
    unittest.main()
