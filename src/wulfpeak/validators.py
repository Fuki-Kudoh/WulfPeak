"""Output-manifest validation that never requires raw inputs or reference files."""

from __future__ import annotations

import json
from pathlib import Path

from .atomic import atomic_write_json
from .models import ValidationResult


class OutputValidationError(ValueError):
    """Raised when validation metadata is missing or malformed."""


BIGWIG_MAGIC = {b"\x26\xfc\x8f\x88", b"\x88\x8f\xfc\x26"}


def validate_path(path: Path, kind: str) -> ValidationResult:
    checks: dict[str, object] = {
        "exists": path.exists(),
        "regular_file": path.is_file(),
    }
    if not path.is_file():
        return ValidationResult(False, checks)
    size = path.stat().st_size
    checks["size"] = size
    warnings: list[str] = []
    if kind in {"narrowPeak", "broadPeak", "bed3"} and size == 0:
        warnings.append("ZERO_PEAKS")
        return ValidationResult(True, checks, tuple(warnings))
    if size == 0:
        return ValidationResult(False, checks)
    if kind == "bigwig":
        with path.open("rb") as handle:
            magic = handle.read(4)
        checks["magic"] = magic.hex()
        return ValidationResult(magic in BIGWIG_MAGIC, checks)
    return ValidationResult(True, checks)


def validate_output_manifest(output_dir: str | Path) -> tuple[Path, bool, dict[str, object]]:
    output = Path(output_dir).expanduser().resolve()
    config_manifest = output / "config" / "manifest.json"
    output_manifest = output / "output_manifest.json"
    if not config_manifest.is_file():
        raise OutputValidationError(f"Missing config manifest: {config_manifest}")
    if not output_manifest.is_file():
        raise OutputValidationError(f"Missing output manifest: {output_manifest}")
    try:
        json.loads(config_manifest.read_text(encoding="utf-8"))
        payload = json.loads(output_manifest.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        raise OutputValidationError(f"Cannot read validation metadata: {exc}") from exc
    entries = payload.get("outputs")
    if not isinstance(entries, list):
        raise OutputValidationError("output_manifest.json must contain an outputs list")
    all_valid = True
    for entry in entries:
        if not isinstance(entry, dict) or "path" not in entry:
            raise OutputValidationError("output manifest contains a malformed entry")
        if not entry.get("expected", True):
            entry.update({"exists": False, "validated": True, "checks": {}, "warnings": []})
            continue
        result = validate_path(Path(str(entry["path"])), str(entry.get("kind", "file")))
        entry.update(
            {
                "exists": bool(result.checks.get("exists")),
                "validated": result.valid,
                "checks": result.checks,
                "warnings": list(result.warnings),
            }
        )
        if entry.get("required", False) and not result.valid:
            all_valid = False
    payload["validated"] = all_valid
    atomic_write_json(output_manifest, payload)
    return output_manifest, all_valid, payload
