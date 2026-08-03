"""Normalized configuration for a WulfPeak run."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

from .models import Assay, ReadLayout


PHASES = (
    "fastqc_raw",
    "trim",
    "fastqc_trimmed",
    "align",
    "bam_process",
    "coverage",
    "peak",
    "pool_bam",
    "pooled_peak",
    "consensus_peak",
    "report",
)

IMPLEMENTED_PHASES = PHASES[: PHASES.index("coverage") + 1]


class ConfigurationError(ValueError):
    """Raised for invalid or unsupported run configuration."""


@dataclass(frozen=True)
class RunConfig:
    samplesheet: Path
    fastq_dir: Path
    output_dir: Path
    assay: Assay
    genome_id: str
    bowtie2_index: Path
    effective_genome_size: int
    threads: int
    read_layout: ReadLayout = ReadLayout.PAIRED_END
    min_mapq: int = 30
    duplicate_policy: str = "remove"
    allow_dovetail: bool = False
    normalization: str = "CPM"
    bin_size: int = 10
    blacklist: Path | None = None
    broad_cutoff: float = 0.1
    keep_intermediates: bool = False
    dry_run: bool = False
    resume: bool = True
    force_from: str | None = None
    stop_after: str | None = None

    def __post_init__(self) -> None:
        errors: list[str] = []
        if not self.genome_id.strip():
            errors.append("--genome-id must be non-empty")
        if self.threads <= 0:
            errors.append("--threads must be a positive integer")
        if self.effective_genome_size <= 0:
            errors.append("--effective-genome-size must be a positive integer")
        if self.min_mapq < 0:
            errors.append("--min-mapq must be zero or greater")
        if self.bin_size <= 0:
            errors.append("--bin-size must be a positive integer")
        if self.duplicate_policy not in {"remove", "keep"}:
            errors.append("--duplicate-policy must be remove or keep")
        if self.normalization not in {"CPM", "RPGC", "None"}:
            errors.append("--normalization must be CPM, RPGC, or None")
        if not 0 < self.broad_cutoff <= 1:
            errors.append("--broad-cutoff must satisfy 0 < value <= 1")
        if self.allow_dovetail and self.read_layout is ReadLayout.SINGLE_END:
            errors.append("--allow-dovetail is only applicable to paired-end runs")
        if self.force_from is not None and self.force_from not in PHASES:
            errors.append(f"unknown --force-from phase: {self.force_from}")
        if self.stop_after is not None and self.stop_after not in PHASES:
            errors.append(f"unknown --stop-after phase: {self.stop_after}")
        if (
            not self.dry_run
            and self.stop_after is not None
            and self.stop_after != "coverage"
        ):
            if self.stop_after in IMPLEMENTED_PHASES:
                errors.append(
                    "non-dry-run execution in this release must stop after coverage"
                )
            else:
                errors.append(
                    f"--stop-after {self.stop_after} is beyond the implemented "
                    "execution boundary (coverage)"
                )
        if (
            not self.dry_run
            and self.force_from is not None
            and self.force_from not in IMPLEMENTED_PHASES
        ):
            errors.append(
                f"--force-from {self.force_from} is beyond the implemented "
                "execution boundary (coverage)"
            )
        if errors:
            raise ConfigurationError("Run configuration failed:\n" + "\n".join(
                f"  - {error}" for error in errors
            ))

    def normalized_options(self) -> dict[str, object]:
        payload = asdict(self)
        payload["samplesheet"] = str(self.samplesheet)
        payload["fastq_dir"] = str(self.fastq_dir)
        payload["output_dir"] = str(self.output_dir)
        payload["bowtie2_index"] = str(self.bowtie2_index)
        payload["blacklist"] = str(self.blacklist) if self.blacklist else None
        payload["assay"] = self.assay.value
        payload["read_layout"] = self.read_layout.value
        return payload
