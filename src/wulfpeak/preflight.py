"""Full-run validation performed before any analysis command starts."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from .config import RunConfig
from .fastq import FastqDiscoveryError, discover_fastq_input
from .manifest import ResolvedSample, build_groups, resolve_sample
from .models import ReplicateGroup, ToolInfo
from .paths import create_output_directories
from .samplesheet import SamplesheetValidationError, load_samplesheet
from .tools import REQUIRED_TOOLS, ToolDetectionError, detect_tools


class PreflightError(ValueError):
    """Raised when one or more complete-run checks fail."""


@dataclass(frozen=True)
class PreflightResult:
    resolved_samples: tuple[ResolvedSample, ...]
    groups: tuple[ReplicateGroup, ...]
    index_files: tuple[Path, ...]
    tools: dict[str, ToolInfo]
    warnings: tuple[dict[str, str], ...]


def validate_bowtie2_index(prefix: str | Path) -> list[Path]:
    base = Path(prefix).expanduser().resolve()
    suffixes = (".1", ".2", ".3", ".4", ".rev.1", ".rev.2")
    families: dict[str, list[Path]] = {
        extension: [Path(f"{base}{suffix}{extension}") for suffix in suffixes]
        for extension in (".bt2", ".bt2l")
    }
    complete = [
        files
        for files in families.values()
        if all(path.is_file() and path.stat().st_size > 0 for path in files)
    ]
    if len(complete) == 1:
        other_extension = ".bt2l" if complete[0][0].name.endswith(".bt2") else ".bt2"
        if any(path.exists() for path in families[other_extension]):
            raise PreflightError("Bowtie2 index must not mix .bt2 and .bt2l families")
        return complete[0]
    existing = [path for files in families.values() for path in files if path.exists()]
    invalid = [
        path
        for files in families.values()
        for path in files
        if not path.is_file() or path.stat().st_size == 0
    ]
    if not existing:
        raise PreflightError(f"No Bowtie2 index files found for prefix {base}")
    raise PreflightError(
        "Incomplete Bowtie2 index; missing or empty files include:\n"
        + "\n".join(f"  - {path}" for path in invalid)
    )


def validate_blacklist(path: Path) -> None:
    if not path.is_file() or not os.access(path, os.R_OK):
        raise PreflightError(f"Blacklist is not a readable regular file: {path}")
    with path.open(encoding="utf-8") as handle:
        for row_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise PreflightError(f"Blacklist row {row_number} has fewer than 3 columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise PreflightError(
                    f"Blacklist row {row_number} has non-integer coordinates"
                ) from exc
            if start < 0 or end < 0 or start >= end:
                raise PreflightError(f"Blacklist row {row_number} has invalid coordinates")


def run_preflight(config: RunConfig) -> PreflightResult:
    try:
        samples = load_samplesheet(config.samplesheet)
    except SamplesheetValidationError as exc:
        raise PreflightError(str(exc)) from exc
    errors: list[str] = []
    resolved: list[ResolvedSample] = []
    for sample in samples:
        try:
            raw = discover_fastq_input(
                sample.input_id, config.fastq_dir, config.read_layout
            )
            for path in (raw.r1, raw.r2):
                if path is None:
                    continue
                if (
                    not path.is_file()
                    or not os.access(path, os.R_OK)
                    or path.stat().st_size == 0
                ):
                    errors.append(
                        f"sample {sample.sample_id!r}: FASTQ must be readable "
                        f"and non-empty: {path}"
                    )
            resolved.append(
                resolve_sample(sample, raw, config.output_dir, config.read_layout)
            )
        except FastqDiscoveryError as exc:
            errors.append(f"sample {sample.sample_id!r}: {exc}")
    if errors:
        raise PreflightError("FASTQ validation failed:\n" + "\n".join(
            f"  - {error}" for error in errors
        ))

    output_paths = [
        path
        for sample in resolved
        for path in (sample.trimmed_fastq.r1, sample.trimmed_fastq.r2)
        if path is not None
    ]
    if len(output_paths) != len(set(output_paths)):
        raise PreflightError("Canonical output path collision detected")

    create_output_directories(config.output_dir)
    if not os.access(config.output_dir, os.W_OK):
        raise PreflightError(f"Output directory is not writable: {config.output_dir}")
    index_files = validate_bowtie2_index(config.bowtie2_index)
    if config.blacklist is not None:
        validate_blacklist(config.blacklist)
    groups = build_groups(samples)
    if (
        config.force_from in {"peak", "pool_bam", "pooled_peak", "consensus_peak"}
        and not groups
    ):
        raise PreflightError(
            f"--force-from {config.force_from} is not applicable because the "
            "samplesheet contains no treatment groups"
        )
    required_tools = tuple(
        name
        for name in REQUIRED_TOOLS
        if groups or name not in {"macs3", "bedtools"}
    )
    try:
        tools, tool_warnings = detect_tools(required_tools)
    except ToolDetectionError as exc:
        raise PreflightError(str(exc)) from exc
    warnings = list(tool_warnings)
    for group in groups:
        if len(group.sample_ids) == 1:
            warnings.append(
                {
                    "code": "SINGLE_REPLICATE_GROUP",
                    "message": f"Group {group.group_id} contains one treatment replicate",
                }
            )
    return PreflightResult(
        tuple(resolved), tuple(groups), tuple(index_files), tools, tuple(warnings)
    )
