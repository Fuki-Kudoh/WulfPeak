from __future__ import annotations

import contextlib
import io
import json
import os
import tempfile
import unittest
from pathlib import Path

from wulfpeak.cli import main
from wulfpeak.tools import REQUIRED_TOOLS


INDEX_SUFFIXES = (".1", ".2", ".3", ".4", ".rev.1", ".rev.2")


@contextlib.contextmanager
def environment(name: str, value: str):
    previous = os.environ.get(name)
    os.environ[name] = value
    try:
        yield
    finally:
        if previous is None:
            os.environ.pop(name, None)
        else:
            os.environ[name] = previous


def make_fake_tools(directory: Path) -> None:
    directory.mkdir()
    for name in REQUIRED_TOOLS:
        executable = directory / name
        executable.write_text(
            f"#!/bin/sh\nprintf '{name} fake-1.0\\n'\n",
            encoding="utf-8",
        )
        executable.chmod(0o755)


def make_index(prefix: Path) -> None:
    for suffix in INDEX_SUFFIXES:
        Path(f"{prefix}{suffix}.bt2").write_bytes(b"index")


def run_args(root: Path, *, single_end: bool) -> list[str]:
    args = [
        "run",
        "--samplesheet",
        str(root / "samples.tsv"),
        "--fastq-dir",
        str(root / "fastq"),
        "--output-dir",
        str(root / "out"),
        "--assay",
        "chipseq",
        "--genome-id",
        "synthetic_reference",
        "--bowtie2-index",
        str(root / "index" / "reference"),
        "--effective-genome-size",
        "1000000",
        "--threads",
        "4",
        "--dry-run",
    ]
    if single_end:
        args.append("--SE")
    return args


class RunCliTests(unittest.TestCase):
    def prepare(self, root: Path) -> Path:
        fake_bin = root / "fake_bin"
        make_fake_tools(fake_bin)
        (root / "fastq").mkdir()
        (root / "index").mkdir()
        make_index(root / "index" / "reference")
        return fake_bin

    def test_single_end_check_writes_null_r2_and_empty_qc_columns(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq = root / "fastq"
            fastq.mkdir()
            read = fastq / "raw_treat_01.fastq.gz"
            read.write_bytes(b"synthetic")
            sheet = root / "samples.tsv"
            sheet.write_text(
                "input_id\tsample_id\tgroup_id\tpeak_type\tcontrol\tqvalue\n"
                "raw_treat_01\ttreat_rep1\tfactor_A\tnarrow\t.\t0.01\n",
                encoding="utf-8",
            )
            self.assertEqual(
                main(
                    [
                        "check",
                        "--samplesheet",
                        str(sheet),
                        "--fastq-dir",
                        str(fastq),
                        "--output-dir",
                        str(root / "out"),
                        "--single-end",
                    ]
                ),
                0,
            )
            manifest = json.loads((root / "out" / "config" / "manifest.json").read_text())
            self.assertEqual(manifest["read_layout"], "single_end")
            self.assertIsNone(manifest["samples"][0]["raw_fastq"]["r2"])
            qc_fields = (
                (root / "out" / "qc" / "sample_qc.tsv")
                .read_text()
                .splitlines()[1]
                .split("\t")
            )
            self.assertEqual(qc_fields[7], "")
            self.assertEqual(qc_fields[9], "")

    def test_single_end_dry_run_writes_plan_without_analysis_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fake_bin = self.prepare(root)
            (root / "samples.tsv").write_text(
                "input_id\tsample_id\tgroup_id\tpeak_type\tcontrol\tqvalue\n"
                "raw_treat_01\ttreat_rep1\tfactor_A\tnarrow\t.\t0.01\n",
                encoding="utf-8",
            )
            (root / "fastq" / "raw_treat_01.fastq.gz").write_bytes(b"synthetic")
            stdout = io.StringIO()
            with environment(
                "PATH", f"{fake_bin}{os.pathsep}{os.environ.get('PATH', '')}"
            ), contextlib.redirect_stdout(stdout):
                code = main(run_args(root, single_end=True))
            self.assertEqual(code, 0, stdout.getvalue())
            plan = json.loads((root / "out" / "config" / "command_plan.json").read_text())
            align = next(step for step in plan["steps"] if step["phase"] == "align")
            self.assertIn("-U", align["commands"][0])
            self.assertNotIn("-1", align["commands"][0])
            trim = next(step for step in plan["steps"] if step["phase"] == "trim")
            self.assertNotIn("--paired", trim["commands"][0])
            self.assertFalse((root / "out" / "bam" / "treat_rep1.final.bam").exists())
            self.assertEqual(
                (root / "out" / "status" / "pipeline.status").read_text().strip(),
                "null",
            )

    def test_paired_end_dry_run_plans_control_and_replicates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fake_bin = self.prepare(root)
            (root / "samples.tsv").write_text(
                "input_id\tsample_id\tgroup_id\tpeak_type\tcontrol\tqvalue\n"
                "raw_treat_01\ttreat_rep1\tfactor_A\tnarrow\tinput_A\t0.01\n"
                "raw_treat_02\ttreat_rep2\tfactor_A\tnarrow\tinput_A\t0.01\n"
                "raw_input_01\tinput_A\tinput_A\tinput\t.\t.\n",
                encoding="utf-8",
            )
            for input_id in ("raw_treat_01", "raw_treat_02", "raw_input_01"):
                (root / "fastq" / f"{input_id}_R1.fastq.gz").write_bytes(b"synthetic-r1")
                (root / "fastq" / f"{input_id}_R2.fastq.gz").write_bytes(b"synthetic-r2")
            with environment("PATH", f"{fake_bin}{os.pathsep}{os.environ.get('PATH', '')}"):
                self.assertEqual(main(run_args(root, single_end=False)), 0)
            plan = json.loads((root / "out" / "config" / "command_plan.json").read_text())
            align = next(step for step in plan["steps"] if step["phase"] == "align")
            self.assertIn("-1", align["commands"][0])
            self.assertIn("-2", align["commands"][0])
            self.assertIn("--no-mixed", align["commands"][0])
            peak = next(step for step in plan["steps"] if step["phase"] == "peak")
            self.assertIn("BAMPE", peak["commands"][0])
            self.assertIn("-c", peak["commands"][0])
            self.assertEqual(
                len(
                    [
                        step
                        for step in plan["steps"]
                        if step["phase"] == "pool_bam"
                    ]
                ),
                1,
            )
            self.assertNotIn("--dovetail", align["commands"][0])
            manifest = json.loads((root / "out" / "config" / "manifest.json").read_text())
            self.assertIn("pooled_bam", manifest["groups"][0]["outputs"])
            metadata = json.loads((root / "out" / "metadata" / "run_metadata.json").read_text())
            self.assertEqual(set(metadata["tools"]), set(REQUIRED_TOOLS))

            dovetail_args = run_args(root, single_end=False)
            dovetail_args.append("--allow-dovetail")
            with environment("PATH", f"{fake_bin}{os.pathsep}{os.environ.get('PATH', '')}"):
                self.assertEqual(main(dovetail_args), 0)
            plan = json.loads((root / "out" / "config" / "command_plan.json").read_text())
            align = next(step for step in plan["steps"] if step["phase"] == "align")
            self.assertIn("--dovetail", align["commands"][0])

    def test_non_dry_run_is_rejected_without_claiming_success(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = run_args(root, single_end=True)
            args.remove("--dry-run")
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                code = main(args)
            self.assertEqual(code, 2)
            self.assertIn("not enabled", stderr.getvalue())
            self.assertFalse((root / "out" / "status" / "pipeline.status").exists())

    def test_dry_run_rejects_empty_fastq_before_tool_detection(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "fastq").mkdir()
            (root / "index").mkdir()
            make_index(root / "index" / "reference")
            (root / "samples.tsv").write_text(
                "input_id\tsample_id\tgroup_id\tpeak_type\tcontrol\tqvalue\n"
                "raw_treat_01\ttreat_rep1\tfactor_A\tnarrow\t.\t0.01\n",
                encoding="utf-8",
            )
            (root / "fastq" / "raw_treat_01.fastq.gz").touch()
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                code = main(run_args(root, single_end=True))
            self.assertEqual(code, 2)
            self.assertIn("non-empty", stderr.getvalue())
            self.assertFalse((root / "out" / "config" / "command_plan.json").exists())

    def test_status_json_does_not_require_inputs_or_tools(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "out"
            (output / "status").mkdir(parents=True)
            (output / "status" / "pipeline.status").write_text("failed\n")
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                self.assertEqual(main(["status", "--output-dir", str(output), "--json"]), 0)
            self.assertEqual(json.loads(stdout.getvalue())["pipeline"], "failed")


if __name__ == "__main__":
    unittest.main()
