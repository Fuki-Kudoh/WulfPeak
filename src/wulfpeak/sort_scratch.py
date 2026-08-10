"""Deterministic, WulfPeak-owned scratch paths for ``samtools sort``."""

from __future__ import annotations

import hashlib
import shutil
import tempfile
from pathlib import Path


def sort_run_namespace(sort_temp_dir: Path, output_dir: Path) -> Path:
    """Return the stable scratch namespace owned by one output directory."""

    canonical_output = str(output_dir.expanduser().resolve())
    digest = hashlib.sha256(canonical_output.encode("utf-8")).hexdigest()[:16]
    return sort_temp_dir.expanduser().resolve() / f"wulfpeak-{digest}"


def sort_sample_directory(
    sort_temp_dir: Path, output_dir: Path, sample_id: str
) -> Path:
    return sort_run_namespace(sort_temp_dir, output_dir) / sample_id


def sort_temp_prefix(sort_temp_dir: Path, output_dir: Path, sample_id: str) -> Path:
    return sort_sample_directory(sort_temp_dir, output_dir, sample_id) / "sort"


def validate_sort_scratch(sort_temp_dir: Path, output_dir: Path) -> Path:
    """Create and genuinely write-test the run-owned scratch namespace."""

    namespace = sort_run_namespace(sort_temp_dir, output_dir)
    probe: Path | None = None
    try:
        namespace.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            dir=namespace,
            prefix=".wulfpeak-write-test-",
            delete=False,
        ) as handle:
            handle.write(b"wulfpeak scratch write test\n")
            handle.flush()
            probe = Path(handle.name)
        probe.unlink()
    except OSError as exc:
        if probe is not None:
            try:
                probe.unlink(missing_ok=True)
            except OSError:
                pass
        raise OSError(
            f"samtools sort scratch is not writable: {sort_temp_dir}: {exc}"
        ) from exc
    return namespace


def prepare_sample_sort_scratch(sample_directory: Path) -> None:
    """Recreate an empty sample-owned namespace before a sort attempt."""

    if sample_directory.exists():
        shutil.rmtree(sample_directory)
    sample_directory.mkdir(parents=True, exist_ok=True)


def cleanup_sample_sort_scratch(sample_directory: Path) -> None:
    """Best-effort removal limited to a WulfPeak-owned sample namespace."""

    try:
        shutil.rmtree(sample_directory)
    except FileNotFoundError:
        pass
    except OSError:
        return
    try:
        sample_directory.parent.rmdir()
    except OSError:
        pass
