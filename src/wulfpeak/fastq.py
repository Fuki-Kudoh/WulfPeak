"""Strict, non-recursive FASTQ discovery for run-level SE or PE layouts."""

from __future__ import annotations

from pathlib import Path

from .models import FastqInput, ReadLayout


class FastqDiscoveryError(ValueError):
    """Base class for FASTQ discovery failures."""


class FastqNotFoundError(FastqDiscoveryError):
    """Raised when no valid input can be found."""


class FastqPairMismatchError(FastqDiscoveryError):
    """Raised when only one mate of a PE input can be found."""


class AmbiguousFastqError(FastqDiscoveryError):
    """Raised when more than one supported input candidate exists."""


# Compatibility name for callers written against v0.1.0.
FastqPair = FastqInput

_EXTENSIONS = ("fastq.gz", "fq.gz")
_R1_STEMS = ("_R1", "_R1_001", "_1")
_R2_STEMS = ("_R2", "_R2_001", "_2")


def _direct_files(fastq_dir: Path) -> dict[str, Path]:
    return {path.name: path.resolve() for path in fastq_dir.iterdir() if path.is_file()}


def _names(input_id: str, stems: tuple[str, ...]) -> list[str]:
    return [
        f"{input_id}{stem}.{extension}"
        for stem in stems
        for extension in _EXTENSIONS
    ]


def _candidates(files: dict[str, Path], names: list[str]) -> list[Path]:
    return [files[name] for name in names if name in files]


def _format_candidates(label: str, candidates: list[Path]) -> str:
    if not candidates:
        return f"{label}: none"
    return f"{label}:\n" + "\n".join(f"  - {path}" for path in candidates)


def discover_fastq_input(
    input_id: str, fastq_dir: str | Path, layout: ReadLayout
) -> FastqInput:
    """Resolve one input using exact supported filenames for ``layout``."""

    directory = Path(fastq_dir).expanduser()
    if not directory.is_dir():
        raise FastqNotFoundError(
            f"FASTQ directory does not exist or is not a directory: {directory.resolve()}"
        )
    files = _direct_files(directory)
    r2_candidates = _candidates(files, _names(input_id, _R2_STEMS))

    if layout is ReadLayout.SINGLE_END:
        single_names = [f"{input_id}.{extension}" for extension in _EXTENSIONS]
        single_names.extend(_names(input_id, _R1_STEMS))
        r1_candidates = _candidates(files, single_names)
        if r2_candidates:
            raise FastqPairMismatchError(
                f"Single-end input_id {input_id!r} has an R2 candidate; WulfPeak "
                f"will not ignore it.\n{_format_candidates('R2 candidates', r2_candidates)}"
            )
        if len(r1_candidates) > 1:
            raise AmbiguousFastqError(
                f"Multiple single-end FASTQ candidates found for input_id {input_id!r}; "
                f"WulfPeak will not guess.\n{_format_candidates('R1 candidates', r1_candidates)}"
            )
        if not r1_candidates:
            raise FastqNotFoundError(
                f"No single-end FASTQ file found for input_id {input_id!r} in "
                f"{directory.resolve()}"
            )
        return FastqInput(r1=r1_candidates[0], r2=None)

    r1_candidates = _candidates(files, _names(input_id, _R1_STEMS))
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
    return FastqInput(r1=r1_candidates[0], r2=r2_candidates[0])


def discover_fastq_pair(input_id: str, fastq_dir: str | Path) -> FastqInput:
    """Backward-compatible paired-end discovery entry point."""

    return discover_fastq_input(input_id, fastq_dir, ReadLayout.PAIRED_END)
