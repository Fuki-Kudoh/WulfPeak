from __future__ import annotations

import contextlib
import io
import json
import os
import tempfile
import unittest
from pathlib import Path

from wulfpeak import __version__
from wulfpeak.cli import main


@contextlib.contextmanager
def working_directory(path: Path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


class CliTests(unittest.TestCase):
    def test_release_version_is_v023(self) -> None:
        self.assertEqual(__version__, "0.2.3")

    def test_check_resolves_default_fastq_dir_and_writes_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            analysis_dir = root / "analysis"
            fastq_dir = root / "fastq"
            analysis_dir.mkdir()
            fastq_dir.mkdir()
            (analysis_dir / "samples.tsv").write_text(
                "input_id\tsample_id\tgroup_id\tpeak_type\tcontrol\tqvalue\n"
                "raw.S1\tcanonical_1\tgroup_1\tbroad\t.\t0.05\n",
                encoding="utf-8",
            )
            r1 = fastq_dir / "raw.S1_R1_001.fastq.gz"
            r2 = fastq_dir / "raw.S1_R2_001.fastq.gz"
            r1.touch()
            r2.touch()
            stdout = io.StringIO()

            with working_directory(analysis_dir), contextlib.redirect_stdout(stdout):
                exit_code = main(["check"])

            self.assertEqual(exit_code, 0)
            self.assertIn(f"R1: {r1.resolve()}", stdout.getvalue())
            self.assertIn(f"R2: {r2.resolve()}", stdout.getvalue())

            output_dir = analysis_dir / "WulfPeak_out"
            manifest = json.loads(
                (output_dir / "config" / "manifest.json").read_text(encoding="utf-8")
            )
            sample = manifest["samples"][0]
            self.assertEqual(
                sample["raw_fastq"],
                {"r1": str(r1.resolve()), "r2": str(r2.resolve())},
            )
            self.assertEqual(
                sample["trimmed_fastq"],
                {
                    "r1": str(
                        (output_dir / "trimmed" / "canonical_1_R1_val_1.fq.gz").resolve()
                    ),
                    "r2": str(
                        (output_dir / "trimmed" / "canonical_1_R2_val_2.fq.gz").resolve()
                    ),
                },
            )

            qc = (output_dir / "qc" / "sample_qc.tsv").read_text(encoding="utf-8")
            self.assertIn(str(r1.resolve()), qc)
            self.assertIn(str(r2.resolve()), qc)
            self.assertIn("canonical_1_R1_val_1.fq.gz", qc)
            self.assertIn("canonical_1_R2_val_2.fq.gz", qc)

    def test_check_reports_all_errors_without_writing_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            samplesheet = root / "samples.tsv"
            samplesheet.write_text(
                "input_id\tsample_id\tgroup_id\tpeak_type\tcontrol\tqvalue\n"
                "missing1\tsample1\tgroup\tnarrow\t.\t0.01\n"
                "missing2\tsample2\tgroup\tnarrow\t.\t0.01\n",
                encoding="utf-8",
            )
            fastq_dir = root / "fastq"
            fastq_dir.mkdir()
            output_dir = root / "out"
            stderr = io.StringIO()

            with contextlib.redirect_stderr(stderr):
                exit_code = main(
                    [
                        "check",
                        "--samplesheet",
                        str(samplesheet),
                        "--fastq-dir",
                        str(fastq_dir),
                        "--output-dir",
                        str(output_dir),
                    ]
                )

            self.assertEqual(exit_code, 2)
            self.assertIn("sample 'sample1'", stderr.getvalue())
            self.assertIn("sample 'sample2'", stderr.getvalue())
            self.assertFalse(output_dir.exists())


if __name__ == "__main__":
    unittest.main()
