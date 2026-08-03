"""Small helpers for crash-resistant metadata writes."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any


def atomic_write_text(path: str | Path, content: str) -> Path:
    """Write text through a sibling temporary file and ``os.replace`` it."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent, prefix=f".{destination.name}.", suffix=".tmp"
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return destination


def atomic_write_json(path: str | Path, payload: Any) -> Path:
    return atomic_write_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")
