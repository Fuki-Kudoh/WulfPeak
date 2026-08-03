"""Write resolved run metadata consumed by downstream WulfPeak steps."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from .fastq import FastqPair
from .paths import canonical_trimmed_fastq_pair
from .samplesheet import Sample


@dataclass(frozen=True)
class ResolvedSample:
    sample: Sample
    raw_fastq: FastqPair
    trimmed_fastq: FastqPair


def resolve_sample(
    sample: Sample, raw_fastq: FastqPair, output_dir: str | Path
) -> ResolvedSample:
    return ResolvedSample(
        sample=sample,
        raw_fastq=raw_fastq,
        trimmed_fastq=canonical_trimmed_fastq_pair(sample.sample_id, output_dir),
    )


def _sample_manifest_record(resolved: ResolvedSample) -> dict[str, object]:
    record: dict[str, object] = asdict(resolved.sample)
    record["raw_fastq"] = {
        "r1": str(resolved.raw_fastq.r1),
        "r2": str(resolved.raw_fastq.r2),
    }
    record["trimmed_fastq"] = {
        "r1": str(resolved.trimmed_fastq.r1),
        "r2": str(resolved.trimmed_fastq.r2),
    }
    return record


def write_manifest(
    resolved_samples: list[ResolvedSample],
    fastq_dir: str | Path,
    output_dir: str | Path,
) -> Path:
    output = Path(output_dir).expanduser().resolve()
    manifest_path = output / "config" / "manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "wulfpeak_version": "0.1.0",
        "read_layout": "paired-end",
        "fastq_dir": str(Path(fastq_dir).expanduser().resolve()),
        "output_dir": str(output),
        "samples": [_sample_manifest_record(sample) for sample in resolved_samples],
    }
    manifest_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return manifest_path


def write_sample_qc(
    resolved_samples: list[ResolvedSample], output_dir: str | Path
) -> Path:
    output = Path(output_dir).expanduser().resolve()
    qc_path = output / "qc" / "sample_qc.tsv"
    qc_path.parent.mkdir(parents=True, exist_ok=True)
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
    with qc_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for resolved in resolved_samples:
            row = asdict(resolved.sample)
            row.update(
                {
                    "fastq_r1": str(resolved.raw_fastq.r1),
                    "fastq_r2": str(resolved.raw_fastq.r2),
                    "trimmed_fastq_r1": str(resolved.trimmed_fastq.r1),
                    "trimmed_fastq_r2": str(resolved.trimmed_fastq.r2),
                }
            )
            writer.writerow(row)
    return qc_path
