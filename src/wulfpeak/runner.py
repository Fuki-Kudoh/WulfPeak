"""Execution-foundation orchestration; analysis execution arrives in PR 2/3."""

from __future__ import annotations

import json
import platform
import socket
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .atomic import atomic_write_json
from .config import ConfigurationError, RunConfig
from .manifest import write_run_manifest
from .plan import build_command_plan
from .preflight import PreflightResult, run_preflight
from .status import RunLock, write_pipeline_status


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


def prepare_dry_run(config: RunConfig, argv: list[str]) -> tuple[Path, Path, Path]:
    if not config.dry_run:
        raise ConfigurationError(
            "Analysis execution is intentionally not enabled in the execution-foundation PR; "
            "use --dry-run to validate and write the complete command plan"
        )
    with RunLock(config.output_dir, argv):
        started = datetime.now(timezone.utc).isoformat()
        result = run_preflight(config)
        ensure_manifest_compatible(config, result)
        manifest_path = write_run_manifest(
            config,
            list(result.resolved_samples),
            list(result.groups),
            list(result.index_files),
        )
        plan = build_command_plan(
            config, result.resolved_samples, result.groups, result.tools
        )
        plan_path = atomic_write_json(
            config.output_dir / "config" / "command_plan.json", plan
        )
        finished = datetime.now(timezone.utc).isoformat()
        metadata = {
            "started_at": started,
            "updated_at": finished,
            "finished_at": finished,
            "argv": argv,
            "cwd": str(Path.cwd()),
            "wulfpeak_version": __version__,
            "python_version": platform.python_version(),
            "hostname": socket.gethostname(),
            "options": config.normalized_options(),
            "tools": _tool_metadata(result),
            "command_plan": str(plan_path),
            "dry_run": True,
            "resume": config.resume,
            "warnings": list(result.warnings),
        }
        metadata_path = atomic_write_json(
            config.output_dir / "metadata" / "run_metadata.json", metadata
        )
        state_path = config.output_dir / "status" / "state.json"
        if not state_path.exists():
            atomic_write_json(state_path, {"steps": {}})
        pipeline_status = config.output_dir / "status" / "pipeline.status"
        if not pipeline_status.exists():
            write_pipeline_status(config.output_dir, "null")
        return manifest_path, plan_path, metadata_path
