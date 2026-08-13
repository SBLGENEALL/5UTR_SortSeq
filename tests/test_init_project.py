import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "init_project.py"


class InitProjectTests(unittest.TestCase):
    def test_lowercase_s_token_and_sample_mapping(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            raw = root / "raw_data"
            config = root / "config"
            raw.mkdir()
            samples = [f"UTR_bin{x}_A{x + 1}" for x in range(1, 7)] + ["UTR_unsorted_H2"]
            for sample in samples:
                for read in ("R1", "R2"):
                    (raw / f"{sample}_s1_{read}_001.fastq.gz").write_bytes(b"")
            (raw / "Undetermined_s0_R1_001.fastq.gz").write_bytes(b"")
            (raw / "Undetermined_s0_R2_001.fastq.gz").write_bytes(b"")

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--raw-dir",
                    str(raw),
                    "--config-dir",
                    str(config),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            with (config / "sample_map.csv").open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(rows[0]["sample_id"], "bin1")
            self.assertEqual(rows[0]["ngs_column"], "UTR_bin1_A2")
            self.assertEqual(rows[-1]["ngs_column"], "UTR_unsorted_H2")


if __name__ == "__main__":
    unittest.main()
