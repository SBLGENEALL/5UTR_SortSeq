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

    def test_actual_indexes_are_auto_imported_from_analysis_sample_sheet(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            raw = root / "raw_data"
            config = root / "config"
            sample_sheet = raw / "Analysis" / "1" / "Data" / "260812_sample_sheet.csv"
            sample_sheet.parent.mkdir(parents=True)
            samples = [f"UTR_bin{x}_A{x + 1}" for x in range(1, 7)] + ["UTR_unsorted_H2"]
            i7_values = [
                "ACGTACGA", "ACGTACGC", "ACGTACGG", "ACGTACGT",
                "TGCATGCA", "TGCATGCC", "TGCATGCG",
            ]
            i5_values = [
                "GATCGATA", "GATCGATC", "GATCGATG", "GATCGATT",
                "CTAGCTAA", "CTAGCTAC", "CTAGCTAG",
            ]
            indexes = []
            for sample, i7, i5 in zip(samples, i7_values, i5_values):
                indexes.append((sample, i7, i5))
                for read in ("R1", "R2"):
                    (raw / f"{sample}_s1_{read}_001.fastq.gz").write_bytes(b"")
            (raw / "Undetermined_s0_R1_001.fastq.gz").write_bytes(b"")
            (raw / "Undetermined_s0_R2_001.fastq.gz").write_bytes(b"")

            with sample_sheet.open("w", newline="", encoding="utf-8") as handle:
                handle.write("[Header]\nFileFormatVersion,2\n\n[BCLConvert_Data]\n")
                writer = csv.writer(handle, lineterminator="\n")
                writer.writerow(["Sample_ID", "Index", "Index2"])
                writer.writerows(indexes)

            completed = subprocess.run(
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
            self.assertIn("Imported and validated i7/i5", completed.stdout)
            generated = (config / "SampleSheet.csv").read_text(encoding="utf-8")
            self.assertNotIn("REPLACE_I7", generated)
            self.assertIn(f"{indexes[0][0]},{indexes[0][1]},{indexes[0][2]}", generated)

    def test_sample_sheet_fastq_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            raw = root / "raw_data"
            raw.mkdir()
            config = root / "config"
            samples = [f"UTR_bin{x}_A{x + 1}" for x in range(1, 7)] + ["UTR_unsorted_H2"]
            for sample in samples:
                (raw / f"{sample}_s1_R1_001.fastq.gz").write_bytes(b"")
            sheet = root / "wrong_sample_sheet.csv"
            wrong_i5 = [
                "TGCATGCA", "TGCATGCC", "TGCATGCG", "TGCATGCT",
                "CATGCATA", "CATGCATC", "CATGCATG",
            ]
            with sheet.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["Sample_ID", "index", "index2"])
                for number, i5 in enumerate(wrong_i5, start=1):
                    writer.writerow([f"WRONG_{number}", "ACGTACGT", i5])
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--raw-dir",
                    str(raw),
                    "--config-dir",
                    str(config),
                    "--run-sample-sheet",
                    str(sheet),
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("must match one-to-one", completed.stderr)


if __name__ == "__main__":
    unittest.main()
