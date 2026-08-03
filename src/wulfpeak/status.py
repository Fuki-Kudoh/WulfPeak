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
    if scope not in {"samples", "groups"}:
        raise ValueError("scope must be 'samples' or 'groups'")
    output = Path(output_dir)
    status_path = output / "status" / scope / scope_id / f"{state.key}.status"
    state_path = output / "status" / "state.json"
    try:
        payload = json.loads(state_path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        payload = {"steps": {}}
    step_id = f"{scope[:-1]}:{scope_id}:{state.key}"
    payload.setdefault("steps", {})[step_id] = asdict(state)
    # Both files use the same atomic helper; JSON is replaced before the compact marker.
    atomic_write_json(state_path, payload)
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
    return {
        "output_dir": str(output),
        "pipeline": pipeline,
        "pipeline_display": pipeline_display,
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
