"""Build and atomically write input-resolution and run manifests."""

from __future__ import annotations

import csv
import hashlib
import io
from dataclasses import asdict, dataclass
from pathlib import Path

from . import __version__
from .atomic import atomic_write_json, atomic_write_text
from .config import RunConfig
from .models import FastqInput, ReadLayout, ReplicateGroup
from .paths import canonical_trimmed_fastq
from .samplesheet import Sample


@dataclass(frozen=True)
class ResolvedSample:
    sample: Sample
    raw_fastq: FastqInput
    trimmed_fastq: FastqInput


def resolve_sample(
    sample: Sample,
    raw_fastq: FastqInput,
    output_dir: str | Path,
    layout: ReadLayout | None = None,
) -> ResolvedSample:
    selected_layout = layout or raw_fastq.layout
    if selected_layout is not raw_fastq.layout:
        raise ValueError("resolved FASTQ layout does not match requested layout")
    return ResolvedSample(
        sample=sample,
        raw_fastq=raw_fastq,
        trimmed_fastq=canonical_trimmed_fastq(
            sample.sample_id, output_dir, selected_layout
        ),
    )


def _fastq_record(fastq: FastqInput) -> dict[str, str | None]:
    return {"r1": str(fastq.r1), "r2": str(fastq.r2) if fastq.r2 else None}


def _sample_manifest_record(
    resolved: ResolvedSample, output_dir: Path, include_outputs: bool
) -> dict[str, object]:
    record: dict[str, object] = asdict(resolved.sample)
    record["raw_fastq"] = _fastq_record(resolved.raw_fastq)
    record["trimmed_fastq"] = _fastq_record(resolved.trimmed_fastq)
    if include_outputs:
        sample_id = resolved.sample.sample_id
        extension = (
            "broadPeak" if resolved.sample.peak_type == "broad" else "narrowPeak"
        )
        record["outputs"] = {
            "bam": str(output_dir / "bam" / f"{sample_id}.final.bam"),
            "bai": str(output_dir / "bam" / f"{sample_id}.final.bam.bai"),
            "bigwig": str(output_dir / "bigwig" / f"{sample_id}.normalized.bw"),
            "peak": (
                None
                if resolved.sample.peak_type == "input"
                else str(
                    output_dir
                    / "peaks"
                    / "replicate"
                    / sample_id
                    / f"{sample_id}.peaks.{extension}"
                )
            ),
        }
    return record


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _base_payload(
    resolved_samples: list[ResolvedSample],
    fastq_dir: str | Path,
    output_dir: str | Path,
    layout: ReadLayout,
    samplesheet: str | Path | None,
    include_outputs: bool,
) -> dict[str, object]:
    output = Path(output_dir).expanduser().resolve()
    payload: dict[str, object] = {
        "schema_version": 2,
        "wulfpeak_version": __version__,
        "read_layout": layout.value,
        "fastq_dir": str(Path(fastq_dir).expanduser().resolve()),
        "output_dir": str(output),
        "samples": [
            _sample_manifest_record(sample, output, include_outputs)
            for sample in resolved_samples
        ],
    }
    if samplesheet is not None:
        sheet = Path(samplesheet).expanduser().resolve()
        payload["samplesheet"] = {"path": str(sheet), "sha256": _sha256(sheet)}
    return payload


def write_manifest(
    resolved_samples: list[ResolvedSample],
    fastq_dir: str | Path,
    output_dir: str | Path,
    layout: ReadLayout = ReadLayout.PAIRED_END,
    samplesheet: str | Path | None = None,
) -> Path:
    """Write the layout-aware input-resolution manifest used by ``check``."""

    output = Path(output_dir).expanduser().resolve()
    payload = _base_payload(
        resolved_samples, fastq_dir, output, layout, samplesheet, False
    )
    return atomic_write_json(output / "config" / "manifest.json", payload)


def build_groups(samples: list[Sample]) -> list[ReplicateGroup]:
    grouped: dict[str, list[Sample]] = {}
    for sample in samples:
        if sample.peak_type != "input":
            grouped.setdefault(sample.group_id, []).append(sample)
    return [
        ReplicateGroup(
            group_id=group_id,
            sample_ids=tuple(sample.sample_id for sample in members),
            peak_type=members[0].peak_type,
            control=members[0].control,
            qvalue=members[0].qvalue,
        )
        for group_id, members in grouped.items()
    ]


def write_run_manifest(
    config: RunConfig,
    resolved_samples: list[ResolvedSample],
    groups: list[ReplicateGroup],
    index_files: list[Path],
) -> Path:
    payload = _base_payload(
        resolved_samples,
        config.fastq_dir,
        config.output_dir,
        config.read_layout,
        config.samplesheet,
        True,
    )
    group_records: list[dict[str, object]] = []
    for group in groups:
        extension = "broadPeak" if group.peak_type == "broad" else "narrowPeak"
        record: dict[str, object] = asdict(group)
        record["outputs"] = {
            "pooled_bam": str(
                config.output_dir
                / "bam"
                / "pooled"
                / f"{group.group_id}.pooled.bam"
            ),
            "pooled_bai": str(
                config.output_dir
                / "bam"
                / "pooled"
                / f"{group.group_id}.pooled.bam.bai"
            ),
            "pooled_peak": str(
                config.output_dir
                / "peaks"
                / "pooled"
                / group.group_id
                / f"{group.group_id}.pooled.peaks.{extension}"
            ),
            "consensus_bed": str(
                config.output_dir
                / "peaks"
                / "consensus"
                / group.group_id
                / f"{group.group_id}.consensus.bed"
            ),
            "consensus_support": str(
                config.output_dir
                / "peaks"
                / "consensus"
                / group.group_id
                / f"{group.group_id}.consensus.support.tsv"
            ),
        }
        group_records.append(record)
    payload.update(
        {
            "assay": config.assay.value,
            "genome": {
                "id": config.genome_id,
                "bowtie2_index_prefix": str(config.bowtie2_index),
                "effective_genome_size": config.effective_genome_size,
                "index_files": [str(path) for path in index_files],
            },
            "parameters": {
                key: value
                for key, value in config.normalized_options().items()
                if key
                not in {
                    "samplesheet",
                    "fastq_dir",
                    "output_dir",
                    "assay",
                    "genome_id",
                    "bowtie2_index",
                    "effective_genome_size",
                    "read_layout",
                }
            },
            "groups": group_records,
        }
    )
    return atomic_write_json(config.output_dir / "config" / "manifest.json", payload)


def write_sample_qc(
    resolved_samples: list[ResolvedSample], output_dir: str | Path
) -> Path:
    output = Path(output_dir).expanduser().resolve()
    qc_path = output / "qc" / "sample_qc.tsv"
    fieldnames = (
        "input_id",
        "sample_id",
        "group_id",
        "peak_type",
        "control",
        "qvalue",
        "fastq_r1",
        "fastq_r2",
        "trimmed_fastq_r1",
        "trimmed_fastq_r2",
    )
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for resolved in resolved_samples:
        row = asdict(resolved.sample)
        row.update(
            {
                "fastq_r1": str(resolved.raw_fastq.r1),
                "fastq_r2": str(resolved.raw_fastq.r2) if resolved.raw_fastq.r2 else "",
                "trimmed_fastq_r1": str(resolved.trimmed_fastq.r1),
                "trimmed_fastq_r2": (
                    str(resolved.trimmed_fastq.r2) if resolved.trimmed_fastq.r2 else ""
                ),
            }
        )
        writer.writerow(row)
    return atomic_write_text(qc_path, buffer.getvalue())
