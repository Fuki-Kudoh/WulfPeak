"""Strict, non-recursive paired-end FASTQ discovery."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


class FastqDiscoveryError(ValueError):
    """Base class for FASTQ discovery failures."""


class FastqNotFoundError(FastqDiscoveryError):
    """Raised when neither mate can be found."""


class FastqPairMismatchError(FastqDiscoveryError):
    """Raised when only one mate of a pair can be found."""


class AmbiguousFastqError(FastqDiscoveryError):
    """Raised when more than one file is a candidate for either mate."""


@dataclass(frozen=True)
class FastqPair:
    """Resolved paths for one paired-end sequencing input."""

    r1: Path
    r2: Path


_EXTENSIONS = ("fastq.gz", "fq.gz")
_MATE_STEMS = {
    "R1": ("_R1", "_R1_001", "_1"),
    "R2": ("_R2", "_R2_001", "_2"),
}


def _candidate_paths(input_id: str, fastq_dir: Path, mate: str) -> list[Path]:
    """Return existing direct children matching the supported names for a mate."""

    expected_names = [
        f"{input_id}{stem}.{extension}"
        for stem in _MATE_STEMS[mate]
        for extension in _EXTENSIONS
    ]
    direct_files = {
        path.name: path.resolve() for path in fastq_dir.iterdir() if path.is_file()
    }
    return [direct_files[name] for name in expected_names if name in direct_files]


def _format_candidates(label: str, candidates: list[Path]) -> str:
    if not candidates:
        return f"{label}: none"
    return f"{label}:\n" + "\n".join(f"  - {path}" for path in candidates)


def discover_fastq_pair(input_id: str, fastq_dir: str | Path) -> FastqPair:
    """Discover exactly one R1/R2 pair for ``input_id``.

    Only supported filenames immediately inside ``fastq_dir`` are considered.
    ``input_id`` is interpolated into exact filenames, so regular-expression and
    glob metacharacters have no special meaning.
    """

    directory = Path(fastq_dir).expanduser()
    if not directory.is_dir():
        raise FastqNotFoundError(
            f"FASTQ directory does not exist or is not a directory: {directory.resolve()}"
        )

    r1_candidates = _candidate_paths(input_id, directory, "R1")
    r2_candidates = _candidate_paths(input_id, directory, "R2")

    if len(r1_candidates) > 1 or len(r2_candidates) > 1:
        raise AmbiguousFastqError(
            f"Multiple FASTQ candidates found for input_id {input_id!r}; "
            "WulfPeak will not guess.\n"
            f"{_format_candidates('R1 candidates', r1_candidates)}\n"
            f"{_format_candidates('R2 candidates', r2_candidates)}"
        )

    if not r1_candidates and not r2_candidates:
        raise FastqNotFoundError(
            f"No paired FASTQ files found for input_id {input_id!r} in "
            f"{directory.resolve()}"
        )

    if not r1_candidates or not r2_candidates:
        missing = "R1" if not r1_candidates else "R2"
        found = r2_candidates if missing == "R1" else r1_candidates
        raise FastqPairMismatchError(
            f"Paired-end FASTQ mismatch for input_id {input_id!r}: {missing} is "
            f"missing; found {found[0]}"
        )

    return FastqPair(r1=r1_candidates[0], r2=r2_candidates[0])
