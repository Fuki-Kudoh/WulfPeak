from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

from wulfpeak.command import CommandExecutionError, CommandRunner
from wulfpeak.consensus import build_consensus_plan, strict_majority_threshold
from wulfpeak.preflight import (
    PreflightError,
    validate_blacklist,
    validate_bowtie2_index,
)
from wulfpeak.plan_schema import artifact_contract, promote_file_artifact
from wulfpeak.signatures import canonical_signature
from wulfpeak.status import (
    LockConflictError,
    RunLock,
    StepState,
    can_reuse_step,
    read_status,
    write_pipeline_status,
    write_step_state,
)
from wulfpeak.validators import validate_output_manifest


INDEX_SUFFIXES = (".1", ".2", ".3", ".4", ".rev.1", ".rev.2")


class FoundationTests(unittest.TestCase):
    def test_accepts_complete_small_and_large_bowtie2_index_families(self) -> None:
        for extension in (".bt2", ".bt2l"):
            with self.subTest(extension=extension), tempfile.TemporaryDirectory() as temporary:
                prefix = Path(temporary) / "reference"
                expected = []
                for suffix in INDEX_SUFFIXES:
                    path = Path(f"{prefix}{suffix}{extension}")
                    path.write_bytes(b"index")
                    expected.append(path)
                self.assertEqual(
                    validate_bowtie2_index(prefix),
                    [path.resolve() for path in expected],
                )

    def test_rejects_incomplete_bowtie2_index(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary) / "reference"
            Path(f"{prefix}.1.bt2").write_bytes(b"index")
            with self.assertRaisesRegex(PreflightError, "Incomplete Bowtie2 index"):
                validate_bowtie2_index(prefix)

    def test_blacklist_requires_valid_bed_like_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            blacklist = Path(temporary) / "blacklist.bed"
            blacklist.write_text("# synthetic\nchr1\t0\t10\textra\n", encoding="utf-8")
            validate_blacklist(blacklist)
            blacklist.write_text("chr1\t10\t5\n", encoding="utf-8")
            with self.assertRaisesRegex(PreflightError, "invalid coordinates"):
                validate_blacklist(blacklist)

    def test_command_runner_checks_return_code_and_logs_argv(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            log = Path(temporary) / "command.log"
            runner = CommandRunner()
            result = runner.run([sys.executable, "-c", "print('ok')"], log)
            self.assertEqual(result.returncode, 0)
            self.assertIn(sys.executable, log.read_text(encoding="utf-8"))
            with self.assertRaises(CommandExecutionError):
                runner.run([sys.executable, "-c", "raise SystemExit(7)"], log)
            self.assertIn("exit_code: 7", log.read_text(encoding="utf-8"))

    def test_command_runner_writes_declared_stdout_atomically(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log = root / "stdout.log"
            destination = root / "artifact.txt"
            destination.write_text("previous\n", encoding="utf-8")
            runner = CommandRunner()
            result = runner.run_spec(
                {
                    "type": "command",
                    "argv": [sys.executable, "-c", "print('replacement')"],
                    "stdout_path": str(destination),
                },
                log,
            )
            self.assertEqual(result.stdout_path, destination)
            self.assertEqual(destination.read_text(), "replacement\n")

            with self.assertRaises(CommandExecutionError):
                runner.run_spec(
                    {
                        "type": "command",
                        "argv": [
                            sys.executable,
                            "-c",
                            "print('partial'); raise SystemExit(6)",
                        ],
                        "stdout_path": str(destination),
                    },
                    log,
                )
            self.assertEqual(destination.read_text(), "replacement\n")
            self.assertEqual(list(root.glob(".artifact.txt.*.tmp")), [])

    def test_command_runner_checks_every_pipeline_return_code(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            log = Path(temporary) / "pipeline.log"
            runner = CommandRunner()
            self.assertEqual(
                runner.run_pipeline(
                    [
                        [sys.executable, "-c", "print('synthetic')"],
                        [sys.executable, "-c", "import sys; sys.stdin.read()"],
                    ],
                    log,
                ),
                [0, 0],
            )
            with self.assertRaises(CommandExecutionError):
                runner.run_pipeline(
                    [
                        [sys.executable, "-c", "print('synthetic'); raise SystemExit(9)"],
                        [sys.executable, "-c", "import sys; sys.stdin.read()"],
                    ],
                    log,
                )
            self.assertIn("exit_codes: [9, 0]", log.read_text(encoding="utf-8"))

    def test_signature_is_stable_for_key_order(self) -> None:
        self.assertEqual(
            canonical_signature({"input": "x", "option": 1}),
            canonical_signature({"option": 1, "input": "x"}),
        )
        self.assertNotEqual(
            canonical_signature({"input": "x", "option": 1}),
            canonical_signature({"input": "x", "option": 2}),
        )

    def test_strict_majority_consensus_plan_for_one_to_five_replicates(self) -> None:
        expected_thresholds = {1: 1, 2: 2, 3: 2, 4: 3, 5: 3}
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            for count, expected in expected_thresholds.items():
                with self.subTest(replicates=count):
                    sample_ids = tuple(f"treat_rep{index}" for index in range(1, count + 1))
                    peaks = tuple(
                        output / "peaks" / f"{sample_id}.narrowPeak"
                        for sample_id in sample_ids
                    )
                    actions, artifacts, metadata = build_consensus_plan(
                        bedtools="/fake/bin/bedtools",
                        output_dir=output,
                        group_id="factor_A",
                        sample_ids=sample_ids,
                        replicate_peaks=peaks,
                    )
                    self.assertEqual(strict_majority_threshold(count), expected)
                    self.assertEqual(metadata["required_count"], expected)
                    self.assertEqual(metadata["sample_ids"], list(sample_ids))
                    self.assertEqual(len(metadata["normalization"]), count)
                    multiinter = next(
                        action
                        for action in actions
                        if action["type"] == "command"
                        and "multiinter" in action["argv"]
                    )
                    self.assertTrue(multiinter["stdout_path"].endswith("multiinter.tsv"))
                    filtering = next(
                        action
                        for action in actions
                        if action.get("operation") == "filter_consensus_support"
                    )
                    self.assertEqual(filtering["minimum_support"], expected)
                    self.assertEqual(filtering["sample_ids"], list(sample_ids))
                    self.assertEqual(
                        {artifact["kind"] for artifact in artifacts}, {"bed3", "tsv"}
                    )
                    for artifact in artifacts:
                        self.assertNotEqual(
                            artifact["temporary_path"], artifact["canonical_path"]
                        )
                        self.assertEqual(
                            artifact["promotion"],
                            "validate_then_atomic_file_replace",
                        )

    def test_file_promotion_replaces_outputs_inside_nonempty_directories(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for directory_name, filenames in (
                ("fastqc_raw", ("read_fastqc.html", "read_fastqc.zip")),
                ("fastqc_trimmed", ("trimmed_fastqc.html", "trimmed_fastqc.zip")),
                ("macs3", ("factor_peaks.narrowPeak", "factor_peaks.xls")),
            ):
                with self.subTest(directory=directory_name):
                    canonical_dir = root / "canonical" / directory_name
                    temporary_dir = root / "temporary" / directory_name
                    canonical_dir.mkdir(parents=True)
                    temporary_dir.mkdir(parents=True)
                    sentinel = canonical_dir / "unrelated.txt"
                    sentinel.write_text("keep\n", encoding="utf-8")
                    artifacts = []
                    for filename in filenames:
                        canonical = canonical_dir / filename
                        temporary_file = temporary_dir / filename
                        canonical.write_text("old\n", encoding="utf-8")
                        temporary_file.write_text("new\n", encoding="utf-8")
                        artifacts.append(
                            artifact_contract(
                                "file",
                                temporary_file,
                                canonical,
                                validator="test",
                            )
                        )
                    for artifact in artifacts:
                        promote_file_artifact(artifact, validated=True)
                    self.assertEqual(sentinel.read_text(encoding="utf-8"), "keep\n")
                    for filename in filenames:
                        self.assertEqual(
                            (canonical_dir / filename).read_text(encoding="utf-8"),
                            "new\n",
                        )
                    self.assertEqual(list(temporary_dir.iterdir()), [])

    def test_directory_artifacts_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "directory artifacts"):
            artifact_contract(
                "directory",
                "/tmp/temporary",
                "/tmp/canonical",
                validator="directory",
            )

    def test_resume_requires_done_signature_outputs_and_upstream(self) -> None:
        state = StepState("trim", "done", signature="abc")
        self.assertTrue(can_reuse_step(state, "abc", outputs_valid=True))
        self.assertFalse(can_reuse_step(state, "different", outputs_valid=True))
        self.assertFalse(can_reuse_step(state, "abc", outputs_valid=False))
        self.assertFalse(
            can_reuse_step(state, "abc", outputs_valid=True, upstream_reusable=False)
        )
        self.assertFalse(can_reuse_step(StepState("trim", "running"), "abc", outputs_valid=True))

    def test_status_and_state_updates_are_machine_readable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            write_pipeline_status(output, "running")
            write_step_state(
                output,
                StepState("trim", "done", signature="abc", output_ids=("sample:s1:trim",)),
                scope="samples",
                scope_id="s1",
            )
            payload = read_status(output)
            self.assertEqual(payload["pipeline"], "running")
            self.assertEqual(payload["steps"]["sample:s1:trim"]["status"], "done")
            self.assertEqual(
                (output / "status" / "samples" / "s1" / "trim.status").read_text().strip(),
                "done",
            )

    def test_status_summarizes_current_phase_and_planned_sample_counts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            (output / "config").mkdir()
            (output / "config" / "command_plan.json").write_text(
                json.dumps(
                    {
                        "steps": [
                            {"scope": "sample", "scope_id": sample, "phase": phase}
                            for phase in ("trim", "align")
                            for sample in ("s1", "s2")
                        ]
                    }
                ),
                encoding="utf-8",
            )
            write_pipeline_status(output, "running")
            write_step_state(
                output,
                StepState("trim", "done", signature="s1-trim"),
                scope="samples",
                scope_id="s1",
            )
            write_step_state(
                output,
                StepState(
                    "trim",
                    "running",
                    signature="s2-trim",
                    started_at="2026-08-07T00:00:00+00:00",
                ),
                scope="samples",
                scope_id="s2",
            )

            payload = read_status(output)
            self.assertEqual(payload["phase"], "trim")
            self.assertEqual(
                payload["phase_counts"],
                {
                    "trim": {"completed": 1, "total": 2},
                    "align": {"completed": 0, "total": 2},
                },
            )
            self.assertEqual(
                payload["running"],
                [
                    {
                        "sample_id": "s2",
                        "phase": "trim",
                        "started_at": "2026-08-07T00:00:00+00:00",
                    }
                ],
            )

    def test_done_status_exposes_invalid_outputs_without_mutation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            write_pipeline_status(output, "done")
            write_step_state(
                output,
                StepState("trim", "done", output_ids=("sample:s1:trim",)),
                scope="samples",
                scope_id="s1",
            )
            (output / "output_manifest.json").write_text(
                json.dumps(
                    {
                        "outputs": [
                            {
                                "id": "sample:s1:trim",
                                "expected": True,
                                "validated": False,
                            }
                        ]
                    }
                ),
                encoding="utf-8",
            )
            payload = read_status(output)
            self.assertEqual(payload["pipeline_display"], "done (invalid output)")
            self.assertEqual(
                payload["steps"]["sample:s1:trim"]["display_status"],
                "done (invalid output)",
            )
            self.assertEqual((output / "status" / "pipeline.status").read_text().strip(), "done")

    def test_run_lock_rejects_concurrent_owner_and_allows_reuse(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with RunLock(temporary, ["first"]):
                with self.assertRaises(LockConflictError):
                    with RunLock(temporary, ["second"]):
                        pass
            with RunLock(temporary, ["third"]):
                pass

    def test_output_validation_updates_manifest_and_accepts_empty_peak(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            (output / "config").mkdir()
            (output / "config" / "manifest.json").write_text("{}\n", encoding="utf-8")
            empty_peak = output / "empty.narrowPeak"
            empty_peak.touch()
            missing = output / "missing.bam"
            (output / "output_manifest.json").write_text(
                json.dumps(
                    {
                        "outputs": [
                            {
                                "id": "sample:s1:peak",
                                "path": str(empty_peak),
                                "kind": "narrowPeak",
                                "required": True,
                                "expected": True,
                            },
                            {
                                "id": "sample:s1:bam",
                                "path": str(missing),
                                "kind": "bam",
                                "required": True,
                                "expected": True,
                            },
                        ]
                    }
                ),
                encoding="utf-8",
            )
            path, valid, payload = validate_output_manifest(output)
            self.assertEqual(path, (output / "output_manifest.json").resolve())
            self.assertFalse(valid)
            self.assertTrue(payload["outputs"][0]["validated"])
            self.assertIn("ZERO_PEAKS", payload["outputs"][0]["warnings"])
            self.assertFalse(payload["outputs"][1]["validated"])


if __name__ == "__main__":
    unittest.main()
