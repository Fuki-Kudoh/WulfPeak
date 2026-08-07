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


def command_argvs(step: dict[str, object]) -> list[list[str]]:
    commands: list[list[str]] = []
    for action in step["actions"]:
        if action["type"] == "command":
            commands.append(action["argv"])
        elif action["type"] == "pipeline":
            commands.extend(command["argv"] for command in action["commands"])
    return commands


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
            align_argv = command_argvs(align)[0]
            self.assertIn("-U", align_argv)
            self.assertNotIn("-1", align_argv)
            trim = next(step for step in plan["steps"] if step["phase"] == "trim")
            self.assertNotIn("--paired", command_argvs(trim)[0])
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
            align_argv = command_argvs(align)[0]
            self.assertIn("-1", align_argv)
            self.assertIn("-2", align_argv)
            self.assertIn("--no-mixed", align_argv)
            peak = next(step for step in plan["steps"] if step["phase"] == "peak")
            peak_argv = command_argvs(peak)[0]
            self.assertIn("BAMPE", peak_argv)
            self.assertIn("-c", peak_argv)
            for phase in ("fastqc_raw", "fastqc_trimmed"):
                fastqc_step = next(
                    step for step in plan["steps"] if step["phase"] == phase
                )
                self.assertEqual(len(fastqc_step["artifacts"]), 4)
                self.assertEqual(
                    {artifact["kind"] for artifact in fastqc_step["artifacts"]},
                    {"html", "zip"},
                )
            raw_macs_artifacts = [
                artifact
                for artifact in peak["artifacts"]
                if Path(artifact["canonical_path"]).parent.name == "macs3"
            ]
            self.assertEqual(len(raw_macs_artifacts), 3)
            self.assertNotIn(
                "directory", {artifact["kind"] for artifact in peak["artifacts"]}
            )
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
            self.assertNotIn("--dovetail", align_argv)
            bam_step = next(
                step for step in plan["steps"] if step["phase"] == "bam_process"
            )
            align_commands = command_argvs(align)
            self.assertEqual(align_commands[1][1], "collate")
            self.assertTrue(
                align["artifacts"][0]["canonical_path"].endswith(
                    ".name-collated.bam"
                )
            )
            bam_pipelines = [
                action
                for action in bam_step["actions"]
                if action["type"] == "pipeline"
            ]
            self.assertEqual(
                [
                    [command["argv"][1] for command in action["commands"]]
                    for action in bam_pipelines
                ],
                [["fixmate", "sort"], ["markdup", "view"]],
            )
            serialized_bam_plan = json.dumps(bam_step)
            for forbidden in ("unsorted.bam", "fixmate.bam", "marked.bam"):
                self.assertNotIn(forbidden, serialized_bam_plan)
            for tool_name in ("flagstat", "stats", "idxstats"):
                action = next(
                    action
                    for action in bam_step["actions"]
                    if action["type"] == "command"
                    and tool_name in action["argv"]
                )
                self.assertIsNotNone(action["stdout_path"])
                self.assertEqual(action["stdout_write"], "atomic_replace")
            for phase in ("bam_process", "coverage", "pool_bam"):
                step = next(item for item in plan["steps"] if item["phase"] == phase)
                for artifact in step["artifacts"]:
                    self.assertNotEqual(
                        artifact["temporary_path"], artifact["canonical_path"]
                    )
                    self.assertEqual(
                        artifact["promotion"],
                        "validate_then_atomic_file_replace",
                    )
            consensus = next(
                step for step in plan["steps"] if step["phase"] == "consensus_peak"
            )
            self.assertEqual(consensus["metadata"]["required_count"], 2)
            multiinter = next(
                action
                for action in consensus["actions"]
                if action["type"] == "command" and "multiinter" in action["argv"]
            )
            self.assertIsNotNone(multiinter["stdout_path"])
            manifest = json.loads((root / "out" / "config" / "manifest.json").read_text())
            self.assertIn("pooled_bam", manifest["groups"][0]["outputs"])
            metadata = json.loads((root / "out" / "metadata" / "run_metadata.json").read_text())
            self.assertEqual(set(metadata["tools"]), set(REQUIRED_TOOLS))

            dovetail_args = run_args(root, single_end=False)
            blacklist = root / "blacklist.bed"
            blacklist.write_text("chr1\t0\t10\n", encoding="utf-8")
            dovetail_args.extend(
                ["--allow-dovetail", "--blacklist", str(blacklist)]
            )
            with environment("PATH", f"{fake_bin}{os.pathsep}{os.environ.get('PATH', '')}"):
                self.assertEqual(main(dovetail_args), 0)
            plan = json.loads((root / "out" / "config" / "command_plan.json").read_text())
            align = next(step for step in plan["steps"] if step["phase"] == "align")
            self.assertIn("--dovetail", command_argvs(align)[0])
            peak = next(step for step in plan["steps"] if step["phase"] == "peak")
            intersect = next(
                action
                for action in peak["actions"]
                if action["type"] == "command" and "intersect" in action["argv"]
            )
            self.assertIsNotNone(intersect["stdout_path"])
            self.assertEqual(intersect["stdout_write"], "atomic_replace")

    def test_non_dry_run_rejects_stop_after_beyond_coverage(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = run_args(root, single_end=True)
            args.remove("--dry-run")
            args.extend(["--stop-after", "peak"])
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                code = main(args)
            self.assertEqual(code, 2)
            self.assertIn("beyond the implemented execution boundary", stderr.getvalue())
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

    def test_dry_run_preserves_existing_done_and_failed_state(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fake_bin = self.prepare(root)
            (root / "samples.tsv").write_text(
                "input_id\tsample_id\tgroup_id\tpeak_type\tcontrol\tqvalue\n"
                "raw_treat_01\ttreat_rep1\tfactor_A\tnarrow\t.\t0.01\n",
                encoding="utf-8",
            )
            (root / "fastq" / "raw_treat_01.fastq.gz").write_bytes(b"synthetic")
            selected_path = f"{fake_bin}{os.pathsep}{os.environ.get('PATH', '')}"
            with environment("PATH", selected_path):
                self.assertEqual(main(run_args(root, single_end=True)), 0)
            state_path = root / "out" / "status" / "state.json"
            pipeline_path = root / "out" / "status" / "pipeline.status"
            for status in ("done", "failed"):
                with self.subTest(status=status):
                    state = {
                        "steps": {
                            "sample:treat_rep1:trim": {
                                "key": "trim",
                                "status": status,
                                "signature": f"preserve-{status}",
                            }
                        }
                    }
                    state_text = json.dumps(state, indent=2) + "\n"
                    state_path.write_text(state_text, encoding="utf-8")
                    pipeline_path.write_text(status + "\n", encoding="utf-8")
                    with environment("PATH", selected_path):
                        self.assertEqual(main(run_args(root, single_end=True)), 0)
                    self.assertEqual(state_path.read_text(encoding="utf-8"), state_text)
                    self.assertEqual(
                        pipeline_path.read_text(encoding="utf-8"), status + "\n"
                    )

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
