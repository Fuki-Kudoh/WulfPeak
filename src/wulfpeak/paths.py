"""Canonical WulfPeak output paths."""

from __future__ import annotations

from pathlib import Path

from .fastq import FastqPair


def canonical_trimmed_fastq_pair(
    sample_id: str, output_dir: str | Path = "WulfPeak_out"
) -> FastqPair:
    """Return sample_id-based Trim Galore output paths used downstream."""

    trimmed_dir = Path(output_dir).expanduser().resolve() / "trimmed"
    return FastqPair(
        r1=trimmed_dir / f"{sample_id}_R1_val_1.fq.gz",
        r2=trimmed_dir / f"{sample_id}_R2_val_2.fq.gz",
    )
