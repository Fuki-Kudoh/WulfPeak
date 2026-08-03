"""Preflight, planning, and coverage-bound per-sample execution."""

from __future__ import annotations

import json
import platform
import socket
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .atomic import atomic_write_json
from .command import CommandExecutionError
from .config import PHASES, ConfigurationError, RunConfig
from .executor import StepExecutor, reset_step_temporary
from .manifest import ResolvedSample, write_run_manifest
from .plan import build_command_plan
from .preflight import PreflightResult, run_preflight
from .signatures import canonical_signature, large_file_fingerprint
from .status import (
    RunLock,
    StepState,
    can_reuse_step,
    write_pipeline_status,
    write_step_state,
)
from .validators import validate_artifact


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _normalized_layout(value: object) -> str | None:
    if value == "paired-end":
        return "paired_end"
    if value in {"paired_end", "single_end"}:
        return str(value)
    return None


def ensure_manifest_compatible(config: RunConfig, result: PreflightResult) -> None:
    path = config.output_dir / "config" / "manifest.json"
    if not path.exists():
        return
    try:
        existing = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ConfigurationError(f"Existing manifest cannot be read: {path}: {exc}") from exc
    schema = existing.get("schema_version", 1)
    if not isinstance(schema, int) or schema > 2:
        raise ConfigurationError(f"Unsupported existing manifest schema: {schema!r}")
    layout = _normalized_layout(existing.get("read_layout"))
    if layout is None:
        raise ConfigurationError("Existing manifest has a missing or unknown read_layout")
    if layout != config.read_layout.value:
        raise ConfigurationError("Cannot resume with a different read layout")
    existing_ids = {
        sample.get("sample_id")
        for sample in existing.get("samples", [])
        if isinstance(sample, dict)
    }
    current_ids = {sample.sample.sample_id for sample in result.resolved_samples}
    if existing_ids != current_ids:
        raise ConfigurationError("Cannot resume with a different sample_id set")
    existing_output = existing.get("output_dir")
    if existing_output and Path(str(existing_output)).resolve() != config.output_dir:
        raise ConfigurationError("Existing manifest belongs to a different output directory")


def _tool_metadata(result: PreflightResult) -> dict[str, object]:
    return {
        name: {
            "path": str(info.path),
            "version": info.version,
            "version_command": list(info.version_command),
        }
        for name, info in result.tools.items()
    }


def _metadata(
    config: RunConfig,
    result: PreflightResult,
    argv: list[str],
    plan_path: Path,
    *,
    started: str,
    finished: str | None,
    completed_through: str | None,
) -> dict[str, object]:
    return {
        "started_at": started,
        "updated_at": finished or _now(),
        "finished_at": finished,
        "argv": argv,
        "cwd": str(Path.cwd()),
        "wulfpeak_version": __version__,
        "python_version": platform.python_version(),
        "hostname": socket.gethostname(),
        "options": config.normalized_options(),
        "tools": _tool_metadata(result),
        "command_plan": str(plan_path),
        "dry_run": config.dry_run,
        "resume": config.resume,
        "completed_through": completed_through,
        "warnings": list(result.warnings),
    }


def _read_states(output_dir: Path) -> dict[str, StepState]:
    try:
        payload = json.loads(
            (output_dir / "status" / "state.json").read_text(encoding="utf-8")
        )
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return {}
    raw_states = payload.get("steps", {})
    if not isinstance(raw_states, dict):
        return {}
    states: dict[str, StepState] = {}
    for step_id, value in raw_states.items():
        if not isinstance(step_id, str) or not isinstance(value, dict):
            continue
        try:
            states[step_id] = StepState(
                key=str(value["key"]),
                status=str(value["status"]),
                signature=(
                    str(value["signature"])
                    if value.get("signature") is not None
                    else None
                ),
                started_at=value.get("started_at"),
                finished_at=value.get("finished_at"),
                exit_code=value.get("exit_code"),
                log_path=value.get("log_path"),
                output_ids=tuple(str(item) for item in value.get("output_ids", [])),
                last_error=value.get("last_error"),
            )
        except (KeyError, TypeError, ValueError):
            continue
    return states


def _artifact_ids(step: dict[str, object]) -> tuple[str, ...]:
    artifacts = step.get("artifacts", [])
    return tuple(
        f"{step['scope']}:{step['scope_id']}:{step['phase']}:{index}"
        for index, _artifact in enumerate(artifacts)
    )


def _state_id(step: dict[str, object]) -> str:
    return f"{step['scope']}:{step['scope_id']}:{step['phase']}"


def _sample_inputs(
    phase: str,
    resolved: ResolvedSample,
    config: RunConfig,
    result: PreflightResult,
) -> list[Path]:
    if phase in {"fastqc_raw", "trim"}:
        paths = [resolved.raw_fastq.r1, resolved.raw_fastq.r2]
    elif phase in {"fastqc_trimmed", "align"}:
        paths = [resolved.trimmed_fastq.r1, resolved.trimmed_fastq.r2]
    elif phase == "bam_process":
        paths = [
            config.output_dir
            / "intermediate"
            / "alignment"
            / f"{resolved.sample.sample_id}.unsorted.bam"
        ]
    elif phase == "coverage":
        bam = config.output_dir / "bam" / f"{resolved.sample.sample_id}.final.bam"
        paths = [bam, Path(str(bam) + ".bai")]
    else:
        paths = []
    selected = [path for path in paths if path is not None]
    if phase == "align":
        selected.extend(result.index_files)
    return selected


def _step_signature(
    step: dict[str, object],
    resolved: ResolvedSample,
    config: RunConfig,
    result: PreflightResult,
) -> str:
    inputs = _sample_inputs(str(step["phase"]), resolved, config, result)
    fingerprints = [large_file_fingerprint(path) for path in inputs]
    return canonical_signature(
        {
            "signature_schema": 1,
            "step": step,
            "inputs": fingerprints,
            "tools": _tool_metadata(result),
        }
    )


def _outputs_valid(step: dict[str, object], *, samtools: Path) -> bool:
    artifacts = step.get("artifacts")
    if not isinstance(artifacts, list):
        return False
    return all(
        isinstance(artifact, dict)
        and validate_artifact(artifact, temporary=False, samtools=samtools).valid
        for artifact in artifacts
    )


def _write_output_manifest(
    config: RunConfig,
    plan: dict[str, object],
    *,
    samtools: Path,
    completed_through: str,
) -> Path:
    outputs: list[dict[str, object]] = []
    for step in plan["steps"]:
        if not isinstance(step, dict):
            continue
        for output_id, artifact in zip(_artifact_ids(step), step["artifacts"]):
            result = validate_artifact(artifact, temporary=False, samtools=samtools)
            outputs.append(
                {
                    "id": output_id,
                    "scope": step["scope"],
                    "scope_id": step["scope_id"],
                    "phase": step["phase"],
                    "path": artifact["canonical_path"],
                    "kind": artifact["kind"],
                    "validator": artifact["validator"],
                    "required": True,
                    "expected": True,
                    "exists": bool(result.checks.get("exists")),
                    "validated": result.valid,
                    "checks": result.checks,
                    "warnings": list(result.warnings),
                }
            )
    valid = all(bool(entry["validated"]) for entry in outputs)
    payload = {
        "schema_version": 2,
        "generated_at": _now(),
        "completed_through": completed_through,
        "validated": valid,
        "outputs": outputs,
    }
    path = atomic_write_json(config.output_dir / "output_manifest.json", payload)
    if not valid:
        raise CommandExecutionError(
            f"Canonical output validation failed while writing {path}"
        )
    return path


def prepare_run(config: RunConfig, argv: list[str]) -> tuple[Path, Path, Path]:
    selected_stop = config.stop_after or (None if config.dry_run else "coverage")
    if not config.dry_run and selected_stop != "coverage":
        raise ConfigurationError(
            "non-dry-run execution must resolve to --stop-after coverage"
        )

    with RunLock(config.output_dir, argv):
        started = _now()
        result = run_preflight(config)
        ensure_manifest_compatible(config, result)
        manifest_path = write_run_manifest(
            config,
            list(result.resolved_samples),
            list(result.groups),
            list(result.index_files),
        )
        plan = build_command_plan(
            config,
            result.resolved_samples,
            result.groups,
            result.tools,
            stop_after=selected_stop,
        )
        plan_path = atomic_write_json(
            config.output_dir / "config" / "command_plan.json", plan
        )
        metadata_path = config.output_dir / "metadata" / "run_metadata.json"

        if config.dry_run:
            finished = _now()
            atomic_write_json(
                metadata_path,
                _metadata(
                    config,
                    result,
                    argv,
                    plan_path,
                    started=started,
                    finished=finished,
                    completed_through=None,
                ),
            )
            state_path = config.output_dir / "status" / "state.json"
            if not state_path.exists():
                atomic_write_json(state_path, {"steps": {}})
            pipeline_status = config.output_dir / "status" / "pipeline.status"
            if not pipeline_status.exists():
                write_pipeline_status(config.output_dir, "null")
            return manifest_path, plan_path, metadata_path

        samtools = result.tools["samtools"].path
        executor = StepExecutor(samtools=samtools)
        states = _read_states(config.output_dir)
        resolved_by_id = {
            sample.sample.sample_id: sample for sample in result.resolved_samples
        }
        chain_unchanged = {
            sample.sample.sample_id: True for sample in result.resolved_samples
        }
        force_index = (
            PHASES.index(config.force_from) if config.force_from is not None else None
        )
        atomic_write_json(
            metadata_path,
            _metadata(
                config,
                result,
                argv,
                plan_path,
                started=started,
                finished=None,
                completed_through=None,
            ),
        )
        write_pipeline_status(config.output_dir, "running")
        try:
            for step in plan["steps"]:
                if not isinstance(step, dict) or step.get("scope") != "sample":
                    raise CommandExecutionError(
                        "The coverage execution boundary only supports sample steps"
                    )
                sample_id = str(step["scope_id"])
                phase = str(step["phase"])
                resolved = resolved_by_id[sample_id]
                signature = _step_signature(step, resolved, config, result)
                log_path = (
                    config.output_dir / "logs" / "samples" / sample_id / f"{phase}.log"
                )
                output_ids = _artifact_ids(step)
                forced = force_index is not None and PHASES.index(phase) >= force_index
                reusable = (
                    config.resume
                    and not forced
                    and can_reuse_step(
                        states.get(_state_id(step)),
                        signature,
                        outputs_valid=_outputs_valid(step, samtools=samtools),
                        upstream_reusable=chain_unchanged[sample_id],
                    )
                )
                if reusable:
                    continue

                chain_unchanged[sample_id] = False
                reset_step_temporary(config.output_dir, step)
                step_started = _now()
                running = StepState(
                    phase,
                    "running",
                    signature=signature,
                    started_at=step_started,
                    log_path=str(log_path),
                    output_ids=output_ids,
                )
                write_step_state(
                    config.output_dir, running, scope="samples", scope_id=sample_id
                )
                try:
                    executor.execute(step, log_path)
                except BaseException as exc:
                    failed = StepState(
                        phase,
                        "failed",
                        signature=signature,
                        started_at=step_started,
                        finished_at=_now(),
                        exit_code=(
                            exc.returncode
                            if isinstance(exc, CommandExecutionError)
                            and exc.returncode is not None
                            else 1
                        ),
                        log_path=str(log_path),
                        output_ids=output_ids,
                        last_error=str(exc),
                    )
                    write_step_state(
                        config.output_dir, failed, scope="samples", scope_id=sample_id
                    )
                    raise CommandExecutionError(
                        f"Step {phase} failed for sample {sample_id}: {exc}"
                    ) from exc
                done = StepState(
                    phase,
                    "done",
                    signature=signature,
                    started_at=step_started,
                    finished_at=_now(),
                    exit_code=0,
                    log_path=str(log_path),
                    output_ids=output_ids,
                )
                write_step_state(
                    config.output_dir, done, scope="samples", scope_id=sample_id
                )
                states[_state_id(step)] = done

            _write_output_manifest(
                config, plan, samtools=samtools, completed_through="coverage"
            )
        except BaseException:
            failed_at = _now()
            write_pipeline_status(config.output_dir, "failed")
            atomic_write_json(
                metadata_path,
                _metadata(
                    config,
                    result,
                    argv,
                    plan_path,
                    started=started,
                    finished=failed_at,
                    completed_through=None,
                ),
            )
            raise

        finished = _now()
        write_pipeline_status(config.output_dir, "done")
        atomic_write_json(
            metadata_path,
            _metadata(
                config,
                result,
                argv,
                plan_path,
                started=started,
                finished=finished,
                completed_through="coverage",
            ),
        )
        return manifest_path, plan_path, metadata_path


def prepare_dry_run(config: RunConfig, argv: list[str]) -> tuple[Path, Path, Path]:
    """Backward-compatible entrypoint now used for both dry and real runs."""

    return prepare_run(config, argv)
