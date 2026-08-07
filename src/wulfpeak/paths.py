"""Canonical WulfPeak output paths and directory setup."""

from __future__ import annotations

from pathlib import Path

from .models import FastqInput, ReadLayout


OUTPUT_DIRECTORIES = (
    "config",
    "metadata",
    "logs/samples",
    "logs/groups",
    "logs/pipeline",
    "status/.locks",
    "status/samples",
    "status/groups",
    "status/pipeline",
    "trimmed",
    "intermediate",
    "bam/pooled",
    "bigwig",
    "peaks/replicate",
    "peaks/pooled",
    "peaks/consensus",
    "qc/fastqc/raw",
    "qc/fastqc/trimmed",
    "qc/trimming",
    "qc/samtools",
    "multiqc",
    "report",
)


def canonical_trimmed_fastq(
    sample_id: str,
    output_dir: str | Path = "WulfPeak_out",
    layout: ReadLayout = ReadLayout.PAIRED_END,
) -> FastqInput:
    trimmed_dir = Path(output_dir).expanduser().resolve() / "trimmed"
    if layout is ReadLayout.SINGLE_END:
        return FastqInput(r1=trimmed_dir / f"{sample_id}_trimmed.fq.gz", r2=None)
    return FastqInput(
        r1=trimmed_dir / f"{sample_id}_R1_val_1.fq.gz",
        r2=trimmed_dir / f"{sample_id}_R2_val_2.fq.gz",
    )


def canonical_trimmed_fastq_pair(
    sample_id: str, output_dir: str | Path = "WulfPeak_out"
) -> FastqInput:
    """Backward-compatible paired-end path helper."""

    return canonical_trimmed_fastq(sample_id, output_dir, ReadLayout.PAIRED_END)


def create_output_directories(output_dir: str | Path) -> Path:
    output = Path(output_dir).expanduser().resolve()
    for relative in OUTPUT_DIRECTORIES:
        (output / relative).mkdir(parents=True, exist_ok=True)
    return output
