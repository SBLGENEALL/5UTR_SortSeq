import csv
import gzip
import importlib.util
import json
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "rescue_dual_index.py"
SPEC = importlib.util.spec_from_file_location("rescue_dual_index", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


SAMPLES = [
    ("bin1", "AAAAAAAA", "CCCCCCCC"),
    ("bin2", "AAAATTTT", "CCCCGGGG"),
    ("bin3", "TTTTAAAA", "GGGGCCCC"),
    ("bin4", "TTTTTTTT", "GGGGGGGG"),
    ("bin5", "ACACACAC", "TGTGTGTG"),
    ("bin6", "CACACACA", "GTGTGTGT"),
    ("unsorted", "AGAGAGAG", "CTCTCTCT"),
]


def write_fastq(path: Path, indexes):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="ascii") as handle:
        for number, index in enumerate(indexes, start=1):
            handle.write(f"@INST:1:FC:1:1:{number}:1 1:N:0:{index}\n")
            handle.write("ACGTACGT\n+\nIIIIIIII\n")


class DualIndexRescueTests(unittest.TestCase):
    def test_actual_filename_style_is_parsed_as_expected(self):
        parsed = MODULE.parse_fastq_name(Path("UTR_bin1_A2_s1_R1_001.fastq.gz"))
        self.assertIsNotNone(parsed)
        self.assertEqual(parsed.sample, "UTR_bin1_A2")
        self.assertEqual(parsed.read, "R1")

    def test_unique_nearest_and_ambiguous_rules(self):
        expected = [MODULE.ExpectedIndex(*row) for row in SAMPLES]
        sample, status, total, d7, d5 = MODULE.classify_index(
            ("AAAATAAA", "CCCCCCCC"), expected, 2, 2, 1
        )
        self.assertEqual(sample, "bin1")
        self.assertEqual((status, total, d7, d5), ("rescued_distance_1", 1, 1, 0))
        sample, status, *_ = MODULE.classify_index(
            ("NNNNNNNN", "NNNNNNNN"), expected, 2, 2, 1
        )
        self.assertIsNone(sample)
        self.assertEqual(status, "too_distant")

    def test_i5_reverse_complement(self):
        expected = [MODULE.ExpectedIndex("x", "AAAAAAAA", "ACGTACGA")]
        oriented = MODULE.orient_expected(expected, "reverse_complement")
        self.assertEqual(oriented[0].i5, "TCGTACGT")

    def test_fast_lookup_is_equivalent_to_reference_classifier(self):
        expected = [MODULE.ExpectedIndex(*row) for row in SAMPLES]
        fast = MODULE.FastIndexClassifier(expected, 2, 2, 1)
        self.assertGreater(len(fast.lookup), 1_000)
        for observed, fast_result in fast.lookup.items():
            reference = MODULE.classify_index(observed, expected, 2, 2, 1)
            self.assertEqual(fast_result, reference)
            self.assertEqual(fast.classify(observed), reference)

        generator = random.Random(260812)
        for _ in range(1_000):
            observed = (
                "".join(generator.choice("ACGTN") for _ in range(8)),
                "".join(generator.choice("ACGTN") for _ in range(8)),
            )
            self.assertEqual(
                fast.classify(observed),
                MODULE.classify_index(observed, expected, 2, 2, 1),
            )

    def test_end_to_end_rescue_writes_extra_chunks_and_residual(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            raw = root / "raw"
            raw.mkdir()
            sample_sheet = root / "SampleSheet.csv"
            with sample_sheet.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["Sample_ID", "index", "index2"])
                writer.writerows(SAMPLES)

            # One existing pair for each assigned sample allows sample-name and
            # i5-orientation validation without rewriting the raw data.
            for sample, i7, i5 in SAMPLES:
                write_fastq(raw / f"{sample}_S1_L001_R1_001.fastq.gz", [f"{i7}+{i5}"])
                write_fastq(raw / f"{sample}_S1_L001_R2_001.fastq.gz", [f"{i7}+{i5}"])

            undetermined_indexes = [
                "AAAAAAAA+CCCCCCCC",  # exact bin1
                "AAAATTTA+CCCCGGGG",  # 1 mismatch bin2
                "TTTTAATA+GGGGCCCA",  # 2 total mismatches bin3
                "NNNNNNNN+NNNNNNNN",  # residual
            ]
            write_fastq(raw / "Undetermined_S0_L001_R1_001.fastq.gz", undetermined_indexes)
            write_fastq(raw / "Undetermined_S0_L001_R2_001.fastq.gz", undetermined_indexes)

            outdir = root / "out"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--input-dir",
                    str(raw),
                    "--sample-sheet",
                    str(sample_sheet),
                    "--outdir",
                    str(outdir),
                    "--orientation-scan-reads",
                    "10",
                    "--progress-interval-seconds",
                    "0",
                    "--progress-check-reads",
                    "1",
                    "--compression-backend",
                    "parallel_python",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            manifest = json.loads((outdir / "rescue_manifest.json").read_text())
            self.assertEqual(manifest["rescued_read_pairs_or_reads"], 3)
            self.assertEqual(manifest["status_counts"]["too_distant"], 1)
            self.assertEqual(manifest["sample_rescued_counts"]["bin1"], 1)
            self.assertEqual(manifest["sample_rescued_counts"]["bin2"], 1)
            self.assertEqual(manifest["sample_rescued_counts"]["bin3"], 1)
            self.assertEqual(
                manifest["performance"]["selected_compression_backend"],
                "parallel_python",
            )
            self.assertGreater(
                manifest["performance"]["index_classifier"]["precomputed_lookup_hits"],
                0,
            )
            self.assertTrue((outdir / "fastq" / "bin1_L900_R1_001.fastq.gz").exists())
            self.assertTrue(
                (outdir / "fastq" / "Undetermined_residual_L900_R1_001.fastq.gz").exists()
            )
            self.assertTrue((outdir / "fastq" / "bin1_S1_L001_R1_001.fastq.gz").is_symlink())
            with gzip.open(
                outdir / "fastq" / "bin1_L900_R1_001.fastq.gz",
                "rt",
                encoding="ascii",
            ) as handle:
                self.assertEqual(sum(1 for _ in handle), 4)
            self.assertIn("[rescue]", completed.stderr)
            self.assertIn("100.00%", completed.stderr)
            progress = json.loads((outdir / "rescue_progress.json").read_text())
            self.assertEqual(progress["status"], "completed")
            self.assertAlmostEqual(progress["progress_percent"], 100.0)
            self.assertEqual(progress["processed_read_pairs_or_reads"], 4)


if __name__ == "__main__":
    unittest.main()
