"""Schema-2 action execution and validate-then-promote step handling."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

from .atomic import atomic_write_text
from .command import CommandRunner
from .plan_schema import promote_artifact
from .validators import validate_artifact


class ArtifactValidationError(ValueError):
    """Raised when a declared temporary output fails validation."""


class StepExecutor:
    def __init__(self, *, samtools: str | Path, runner: CommandRunner | None = None):
        self.samtools = str(samtools)
        self.runner = runner or CommandRunner()

    def execute(self, step: dict[str, object], log_path: str | Path) -> None:
        actions = step.get("actions")
        if not isinstance(actions, list):
            raise ValueError("step actions must be a list")
        for action in actions:
            if not isinstance(action, dict):
                raise ValueError("step action must be an object")
            action_type = action.get("type")
            if action_type == "command":
                self.runner.run_spec(action, log_path)
            elif action_type == "pipeline":
                commands = action.get("commands")
                if not isinstance(commands, list):
                    raise ValueError("pipeline commands must be a list")
                argvs: list[list[str]] = []
                for command in commands:
                    if not isinstance(command, dict):
                        raise ValueError("pipeline command must be an object")
                    argv = command.get("argv")
                    if not isinstance(argv, list) or not all(
                        isinstance(item, str) for item in argv
                    ):
                        raise ValueError("pipeline argv must be a list of strings")
                    argvs.append(argv)
                self.runner.run_pipeline(argvs, log_path)
            elif action_type == "internal":
                self._run_internal(action)
            else:
                raise ValueError(f"unsupported action type: {action_type!r}")

        artifacts = step.get("artifacts")
        if not isinstance(artifacts, list) or not all(
            isinstance(artifact, dict) for artifact in artifacts
        ):
            raise ValueError("step artifacts must be a list of objects")
        results = [
            validate_artifact(artifact, temporary=True, samtools=self.samtools)
            for artifact in artifacts
        ]
        failures = [
            f"{artifact.get('temporary_path')}: {result.checks}"
            for artifact, result in zip(artifacts, results)
            if not result.valid
        ]
        if failures:
            raise ArtifactValidationError(
                "Temporary artifact validation failed:\n"
                + "\n".join(f"  - {failure}" for failure in failures)
            )
        for artifact in artifacts:
            promote_artifact(artifact, validated=True)

    def _run_internal(self, action: dict[str, object]) -> None:
        operation = action.get("operation")
        if operation != "normalize_trim_galore_outputs":
            raise ValueError(f"unsupported internal operation: {operation!r}")
        output_dir = Path(str(action["output_dir"]))
        layout = str(action["read_layout"])
        targets_value = action.get("canonical_temporary_paths")
        if not isinstance(targets_value, list) or not all(
            isinstance(item, str) for item in targets_value
        ):
            raise ValueError("canonical_temporary_paths must be a list of strings")
        targets = [Path(item) for item in targets_value]
        fastq_targets = [path for path in targets if path.name.endswith(".fq.gz")]
        report_targets = [path for path in targets if path.name.endswith(".txt")]
        expected_suffixes = (
            ("_val_1.fq.gz", "_val_2.fq.gz")
            if layout == "paired_end"
            else ("_trimmed.fq.gz",)
        )
        if len(fastq_targets) != len(expected_suffixes) or len(report_targets) != 1:
            raise ValueError("Trim Galore normalization targets do not match read layout")

        sources: list[Path] = []
        for suffix in expected_suffixes:
            matches = sorted(
                path
                for path in output_dir.glob(f"*{suffix}")
                if path.is_file() and "canonical" not in path.parts
            )
            if len(matches) != 1:
                raise ValueError(
                    f"Expected exactly one Trim Galore output matching *{suffix}; "
                    f"found {len(matches)}"
                )
            sources.append(matches[0])
        for source, target in zip(sources, fastq_targets):
            target.parent.mkdir(parents=True, exist_ok=True)
            os.replace(source, target)

        reports = sorted(
            path
            for path in output_dir.glob("*_trimming_report.txt")
            if path.is_file() and "canonical" not in path.parts
        )
        if len(reports) != len(expected_suffixes):
            raise ValueError(
                "Expected one Trim Galore report per input read; "
                f"found {len(reports)}"
            )
        report_text = ""
        for report in reports:
            if len(reports) > 1:
                report_text += f"## {report.name}\n"
            report_text += report.read_text(encoding="utf-8")
            if report_text and not report_text.endswith("\n"):
                report_text += "\n"
        atomic_write_text(report_targets[0], report_text)


def reset_step_temporary(output_dir: Path, step: dict[str, object]) -> Path:
    scope = str(step["scope"])
    scope_id = str(step["scope_id"])
    phase = str(step["phase"])
    scope_dir = {"sample": "samples", "group": "groups", "pipeline": "pipeline"}.get(
        scope
    )
    if scope_dir is None:
        raise ValueError(f"unsupported step scope: {scope!r}")
    temporary = output_dir / "intermediate" / ".steps" / scope_dir / scope_id / phase
    if temporary.exists():
        shutil.rmtree(temporary)
    temporary.mkdir(parents=True, exist_ok=True)
    return temporary
