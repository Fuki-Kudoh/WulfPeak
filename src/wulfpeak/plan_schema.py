"""Structured dry-run actions and temporary-to-canonical artifact contracts."""

from __future__ import annotations

from pathlib import Path


def command_action(
    argv: list[str], *, stdout_path: str | Path | None = None
) -> dict[str, object]:
    return {
        "type": "command",
        "argv": argv,
        "stdout_path": str(stdout_path) if stdout_path is not None else None,
        "stdout_write": "atomic_replace" if stdout_path is not None else None,
    }


def pipeline_action(commands: list[list[str]]) -> dict[str, object]:
    return {
        "type": "pipeline",
        "commands": [{"argv": command} for command in commands],
        "check_all_return_codes": True,
    }


def internal_action(operation: str, **parameters: object) -> dict[str, object]:
    return {"type": "internal", "operation": operation, **parameters}


def artifact_contract(
    kind: str,
    temporary_path: str | Path,
    canonical_path: str | Path,
    *,
    validator: str,
) -> dict[str, object]:
    temporary = Path(temporary_path)
    canonical = Path(canonical_path)
    if temporary == canonical:
        raise ValueError("temporary and canonical artifact paths must differ")
    return {
        "kind": kind,
        "temporary_path": str(temporary),
        "canonical_path": str(canonical),
        "validator": validator,
        "promotion": "validate_then_atomic_replace",
    }
