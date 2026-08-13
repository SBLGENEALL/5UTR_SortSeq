import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "run_pipeline.sh"


class PipelineStatusTests(unittest.TestCase):
    def test_resources_reports_parallel_settings(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            config = root / "sortseq.env"
            config.write_text(
                f"RESULTS_DIR={root / 'results'}\n"
                "WORKERS=128\n"
                "LIBRARYQC_BATCH_SIZE=10000\n"
                "RESCUE_COMPRESSION_BACKEND=parallel_python\n",
                encoding="utf-8",
            )
            completed = subprocess.run(
                ["bash", str(SCRIPT), "resources"],
                env={**os.environ, "SORTSEQ_PROJECT_CONFIG": str(config)},
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertIn("libraryqc_workers: 128", completed.stdout)
            self.assertIn("libraryqc_batch_size: 10000", completed.stdout)
            self.assertIn("rescue_compression: parallel_python", completed.stdout)

    def test_status_reports_rescue_progress_and_checkpoints(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            results = root / "results"
            rescue = results / "index_rescue"
            rescue.mkdir(parents=True)
            (rescue / "rescue_inspection.json").write_text("{}\n", encoding="utf-8")
            (rescue / "rescue_progress.json").write_text(
                json.dumps(
                    {
                        "status": "running",
                        "progress_percent": 42.5,
                        "chunk": 1,
                        "chunks_total": 1,
                        "processed_read_pairs_or_reads": 1_250_000,
                        "rescued_read_pairs_or_reads": 1_000_000,
                        "rescued_percent_so_far": 80.0,
                        "read_pairs_or_reads_per_second": 25_000,
                        "eta_seconds": 120,
                        "updated_at": "2026-08-13T00:00:00+00:00",
                    }
                )
                + "\n",
                encoding="utf-8",
            )
            config = root / "sortseq.env"
            config.write_text(
                f"RESULTS_DIR={results}\nPYTHON_BIN={os.environ.get('PYTHON', 'python3')}\n",
                encoding="utf-8",
            )
            completed = subprocess.run(
                ["bash", str(SCRIPT), "status"],
                env={**os.environ, "SORTSEQ_PROJECT_CONFIG": str(config)},
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertIn("progress: 42.50%", completed.stdout)
            self.assertIn("processed: 1,250,000", completed.stdout)
            self.assertIn("[DONE] 1/5 preflight", completed.stdout)
            self.assertIn("[----] 2/5 rescue", completed.stdout)


if __name__ == "__main__":
    unittest.main()
