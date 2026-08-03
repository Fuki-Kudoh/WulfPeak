"""Stable signatures used by later execution steps for conservative resume."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


def canonical_signature(payload: dict[str, object]) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def large_file_fingerprint(path: str | Path) -> dict[str, object]:
    resolved = Path(path).resolve()
    stat = resolved.stat()
    return {"path": str(resolved), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
