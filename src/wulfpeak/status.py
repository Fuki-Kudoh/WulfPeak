"""Atomic status/state helpers and the single-run advisory lock."""

from __future__ import annotations

import fcntl
import json
import os
import socket
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import IO

from .atomic import atomic_write_json, atomic_write_text
from .config import PHASES


STATUS_VALUES = frozenset({"null", "running", "done", "failed"})


class LockConflictError(RuntimeError):
    """Raised when another process owns the output-directory run lock."""


@dataclass(frozen=True)
class StepState:
    key: str
    status: str
    signature: str | None = None
    started_at: str | None = None
    finished_at: str | None = None
    exit_code: int | None = None
    log_path: str | None = None
    output_ids: tuple[str, ...] = ()
    last_error: str | None = None

    def __post_init__(self) -> None:
        if self.status not in STATUS_VALUES:
            raise ValueError(f"invalid status: {self.status}")


def write_pipeline_status(output_dir: str | Path, status: str) -> Path:
    if status not in STATUS_VALUES:
        raise ValueError(f"invalid status: {status}")
    return atomic_write_text(
        Path(output_dir) / "status" / "pipeline.status", status + "\n"
    )


def write_step_state(
    output_dir: str | Path, state: StepState, *, scope: str, scope_id: str
) -> None:
    state_scopes = {"samples": "sample", "groups": "group", "pipeline": "pipeline"}
    if scope not in state_scopes:
        raise ValueError("scope must be 'samples', 'groups', or 'pipeline'")
    output = Path(output_dir)
    status_path = output / "status" / scope / scope_id / f"{state.key}.status"
    state_path = output / "status" / "state.json"
    step_id = f"{state_scopes[scope]}:{scope_id}:{state.key}"
    lock_path = output / "status" / ".locks" / "state.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    # Use a separate open file description per call: flock then coordinates both
    # worker threads today and independent processes in a future scheduler.
    with lock_path.open("a+", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            try:
                payload = json.loads(state_path.read_text(encoding="utf-8"))
            except FileNotFoundError:
                payload = {"steps": {}}
            steps = payload.setdefault("steps", {})
            if not isinstance(steps, dict):
                raise ValueError(f"invalid shared step state in {state_path}")
            steps[step_id] = asdict(state)
            atomic_write_json(state_path, payload)
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
    # Compact markers have unique sample/phase paths and need no shared lock.
    atomic_write_text(status_path, state.status + "\n")


def read_status(output_dir: str | Path) -> dict[str, object]:
    output = Path(output_dir).expanduser().resolve()
    pipeline_path = output / "status" / "pipeline.status"
    try:
        pipeline = pipeline_path.read_text(encoding="utf-8").strip()
    except FileNotFoundError:
        pipeline = "null"
    state_path = output / "status" / "state.json"
    try:
        state = json.loads(state_path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        state = {"steps": {}}
    invalid_ids: set[str] = set()
    try:
        output_payload = json.loads(
            (output / "output_manifest.json").read_text(encoding="utf-8")
        )
        for entry in output_payload.get("outputs", []):
            if (
                isinstance(entry, dict)
                and entry.get("expected", True)
                and not entry.get("validated", False)
            ):
                invalid_ids.add(str(entry.get("id")))
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        pass
    steps = state.get("steps", {})
    if isinstance(steps, dict):
        for step in steps.values():
            if not isinstance(step, dict) or step.get("status") != "done":
                continue
            output_ids = {str(item) for item in step.get("output_ids", [])}
            step["display_status"] = (
                "done (invalid output)" if output_ids & invalid_ids else "done"
            )
    pipeline_display = (
        "done (invalid output)" if pipeline == "done" and invalid_ids else pipeline
    )
    phase_step_ids: dict[str, list[str]] = {}
    try:
        plan_payload = json.loads(
            (output / "config" / "command_plan.json").read_text(encoding="utf-8")
        )
        planned_steps = plan_payload.get("steps", [])
        if isinstance(planned_steps, list):
            for step in planned_steps:
                if not isinstance(step, dict):
                    continue
                scope = step.get("scope")
                phase = step.get("phase")
                scope_id = step.get("scope_id")
                if (
                    scope in {"sample", "group", "pipeline"}
                    and isinstance(phase, str)
                    and isinstance(scope_id, str)
                ):
                    phase_step_ids.setdefault(phase, []).append(
                        f"{scope}:{scope_id}:{phase}"
                    )
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        pass
    if not phase_step_ids and isinstance(steps, dict):
        for step_id, step in steps.items():
            if not isinstance(step_id, str):
                continue
            parts = step_id.split(":", 2)
            if (
                len(parts) != 3
                or parts[0] not in {"sample", "group", "pipeline"}
                or not isinstance(step, dict)
                or not isinstance(step.get("key"), str)
            ):
                continue
            phase_step_ids.setdefault(str(step["key"]), []).append(step_id)

    phase_order = sorted(
        phase_step_ids,
        key=lambda phase: (
            PHASES.index(phase) if phase in PHASES else len(PHASES),
            phase,
        ),
    )
    phase_counts: dict[str, dict[str, int]] = {}
    if isinstance(steps, dict):
        for phase in phase_order:
            step_ids = phase_step_ids[phase]
            phase_counts[phase] = {
                "completed": sum(
                    1
                    for step_id in step_ids
                    if isinstance(steps.get(step_id), dict)
                    and steps[step_id].get("display_status") == "done"
                ),
                "total": len(step_ids),
            }

    running: list[dict[str, object]] = []
    failed: list[dict[str, object]] = []
    if isinstance(steps, dict):
        for step_id, step in steps.items():
            if not isinstance(step_id, str) or not isinstance(step, dict):
                continue
            parts = step_id.split(":", 2)
            if len(parts) != 3 or parts[0] not in {"sample", "group", "pipeline"}:
                continue
            detail: dict[str, object] = {
                "phase": parts[2],
                "started_at": step.get("started_at"),
            }
            if parts[0] == "sample":
                detail["sample_id"] = parts[1]
            else:
                detail.update({"scope": parts[0], "scope_id": parts[1]})
            if step.get("status") == "running":
                running.append(detail)
            elif step.get("status") == "failed":
                failed.append({**detail, "last_error": step.get("last_error")})
    phase_rank = {phase: index for index, phase in enumerate(phase_order)}
    def detail_key(detail: dict[str, object]) -> tuple[int, str]:
        return (
            phase_rank.get(str(detail["phase"]), len(phase_rank)),
            str(detail.get("sample_id", detail.get("scope_id", ""))),
        )

    running.sort(key=detail_key)
    failed.sort(key=detail_key)

    current_phase: str | None = None
    if running:
        current_phase = str(running[0]["phase"])
    elif pipeline == "failed" and failed:
        current_phase = str(failed[0]["phase"])
    elif pipeline in {"running", "failed"}:
        current_phase = next(
            (
                phase
                for phase, counts in phase_counts.items()
                if counts["completed"] < counts["total"]
            ),
            None,
        )
    return {
        "output_dir": str(output),
        "pipeline": pipeline,
        "pipeline_display": pipeline_display,
        "phase": current_phase,
        "phase_counts": phase_counts,
        "running": running,
        "failed": failed,
        **state,
    }


def can_reuse_step(
    state: StepState | None,
    current_signature: str,
    *,
    outputs_valid: bool,
    upstream_reusable: bool = True,
) -> bool:
    """Apply the conservative five-part skip rule without mutating state."""

    return bool(
        state is not None
        and state.status == "done"
        and state.signature == current_signature
        and outputs_valid
        and upstream_reusable
    )


class RunLock:
    """Non-blocking advisory lock whose file may safely remain after release."""

    def __init__(self, output_dir: str | Path, argv: list[str] | None = None):
        self.path = (
            Path(output_dir).expanduser().resolve()
            / "status"
            / ".locks"
            / "run.lock"
        )
        self.argv = list(argv if argv is not None else sys.argv)
        self._handle: IO[str] | None = None

    def __enter__(self) -> "RunLock":
        self.path.parent.mkdir(parents=True, exist_ok=True)
        handle = self.path.open("a+", encoding="utf-8")
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            handle.seek(0)
            owner = handle.read().strip() or "owner metadata unavailable"
            handle.close()
            raise LockConflictError(
                f"Another WulfPeak run holds {self.path}: {owner}"
            ) from exc
        metadata = {
            "pid": os.getpid(),
            "hostname": socket.gethostname(),
            "started_at": datetime.now(timezone.utc).isoformat(),
            "argv": self.argv,
        }
        handle.seek(0)
        handle.truncate()
        json.dump(metadata, handle, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
        self._handle = handle
        return self

    def __exit__(self, exc_type, exc, traceback) -> None:
        if self._handle is not None:
            fcntl.flock(self._handle.fileno(), fcntl.LOCK_UN)
            self._handle.close()
            self._handle = None
