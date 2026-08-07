"""Deterministic command and artifact plans for the execution-foundation PR."""

from __future__ import annotations

from pathlib import Path

from .config import PHASES, RunConfig
from .consensus import build_consensus_plan
from .manifest import ResolvedSample
from .models import ReadLayout, ReplicateGroup, ToolInfo
from .plan_schema import (
    artifact_contract,
    command_action,
    internal_action,
    pipeline_action,
)
from .resources import allocate_threads


def _tool(tools: dict[str, ToolInfo], name: str) -> str:
    return str(tools[name].path)


def _step_tmp(output: Path, scope: str, scope_id: str, phase: str) -> Path:
    return output / "intermediate" / ".steps" / scope / scope_id / phase


def _fastqc_basename(read: str | Path) -> str:
    name = Path(read).name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    raise ValueError(f"unsupported FASTQ filename in command plan: {name}")


def _fastqc_artifacts(
    reads: list[str], temporary_dir: Path, canonical_dir: Path
) -> list[dict[str, object]]:
    artifacts: list[dict[str, object]] = []
    for read in reads:
        basename = _fastqc_basename(read)
        for extension, kind, validator in (
            ("html", "html", "fastqc_html"),
            ("zip", "zip", "fastqc_zip"),
        ):
            filename = f"{basename}_fastqc.{extension}"
            artifacts.append(
                artifact_contract(
                    kind,
                    temporary_dir / filename,
                    canonical_dir / filename,
                    validator=validator,
                )
            )
    return artifacts


def _record(
    phase: str,
    scope: str,
    scope_id: str,
    actions: list[dict[str, object]],
    artifacts: list[dict[str, object]],
    *,
    metadata: dict[str, object] | None = None,
) -> dict[str, object]:
    record: dict[str, object] = {
        "phase": phase,
        "scope": scope,
        "scope_id": scope_id,
        "actions": actions,
        "artifacts": artifacts,
        "outputs": [artifact["canonical_path"] for artifact in artifacts],
    }
    if metadata:
        record["metadata"] = metadata
    return record


def _phase_major(steps: list[dict[str, object]]) -> list[dict[str, object]]:
    """Keep input order within each phase while placing phase barriers in the plan."""

    return sorted(steps, key=lambda step: PHASES.index(str(step["phase"])))


def _peak_plan(
    *,
    config: RunConfig,
    tools: dict[str, ToolInfo],
    scope: str,
    scope_id: str,
    name: str,
    treatment_bam: Path,
    control_bam: Path | None,
    peak_type: str,
    qvalue: str,
    canonical_root: Path,
    canonical_name: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    step_tmp = _step_tmp(config.output_dir, f"{scope}s", scope_id, "peak")
    if scope == "group":
        step_tmp = _step_tmp(
            config.output_dir, "groups", scope_id, "pooled_peak"
        )
    macs_tmp = step_tmp / "macs3"
    extension = "broadPeak" if peak_type == "broad" else "narrowPeak"
    raw_peak = macs_tmp / f"{name}_peaks.{extension}"
    canonical_peak_tmp = step_tmp / "canonical" / canonical_name
    canonical_peak = canonical_root / canonical_name
    macs = [
        _tool(tools, "macs3"),
        "callpeak",
        "-t",
        str(treatment_bam),
    ]
    if control_bam is not None:
        macs.extend(["-c", str(control_bam)])
    macs.extend(
        [
            "-f",
            "BAMPE" if config.read_layout is ReadLayout.PAIRED_END else "BAM",
            "-g",
            str(config.effective_genome_size),
            "-n",
            name,
            "-q",
            qvalue,
            "--outdir",
            str(macs_tmp),
        ]
    )
    if peak_type == "broad":
        macs.extend(["--broad", "--broad-cutoff", str(config.broad_cutoff)])
    actions = [command_action(macs)]
    if config.blacklist is not None:
        actions.append(
            command_action(
                [
                    _tool(tools, "bedtools"),
                    "intersect",
                    "-v",
                    "-a",
                    str(raw_peak),
                    "-b",
                    str(config.blacklist),
                ],
                stdout_path=canonical_peak_tmp,
            )
        )
    else:
        actions.append(
            internal_action(
                "copy_peak_to_canonical_temp",
                source_path=str(raw_peak),
                destination_path=str(canonical_peak_tmp),
                write_mode="atomic_replace",
            )
        )
    raw_output_specs = [
        (f"{name}_peaks.{extension}", extension, f"{extension}_or_empty"),
        (f"{name}_peaks.xls", "text", "text_nonempty"),
    ]
    if peak_type == "broad":
        raw_output_specs.append(
            (f"{name}_peaks.gappedPeak", "gappedPeak", "gappedPeak_or_empty")
        )
    else:
        raw_output_specs.append(
            (f"{name}_summits.bed", "bed", "bed_or_empty")
        )
    artifacts = [
        artifact_contract(
            kind,
            macs_tmp / filename,
            canonical_root / "macs3" / filename,
            validator=validator,
        )
        for filename, kind, validator in raw_output_specs
    ]
    artifacts.append(
        artifact_contract(
            extension,
            canonical_peak_tmp,
            canonical_peak,
            validator=f"{extension}_or_empty",
        )
    )
    return actions, artifacts


def build_command_plan(
    config: RunConfig,
    samples: tuple[ResolvedSample, ...],
    groups: tuple[ReplicateGroup, ...],
    tools: dict[str, ToolInfo],
    *,
    stop_after: str | None = None,
) -> dict[str, object]:
    """Return complete actions and promotion contracts without executing them."""

    output = config.output_dir
    thread_allocation = allocate_threads(config.threads)
    selected_stop = stop_after or PHASES[-1]
    stop_index = PHASES.index(selected_stop)

    def includes(phase: str) -> bool:
        return PHASES.index(phase) <= stop_index

    steps: list[dict[str, object]] = []
    for resolved in samples:
        sample = resolved.sample
        sample_id = sample.sample_id
        raw_reads = [str(resolved.raw_fastq.r1)]
        trimmed_reads = [str(resolved.trimmed_fastq.r1)]
        if resolved.raw_fastq.r2 is not None:
            raw_reads.append(str(resolved.raw_fastq.r2))
        if resolved.trimmed_fastq.r2 is not None:
            trimmed_reads.append(str(resolved.trimmed_fastq.r2))

        raw_qc_tmp = _step_tmp(output, "samples", sample_id, "fastqc_raw")
        raw_qc = output / "qc" / "fastqc" / "raw" / sample_id
        steps.append(
            _record(
                "fastqc_raw",
                "sample",
                sample_id,
                [
                    command_action(
                        [
                            _tool(tools, "fastqc"),
                            "--threads",
                            str(config.threads),
                            "--outdir",
                            str(raw_qc_tmp),
                            *raw_reads,
                        ]
                    )
                ],
                _fastqc_artifacts(raw_reads, raw_qc_tmp, raw_qc),
            )
        )

        trim_tmp = _step_tmp(output, "samples", sample_id, "trim")
        trim_command = [
            _tool(tools, "trim_galore"),
            "--cores",
            str(config.threads),
        ]
        if config.read_layout is ReadLayout.PAIRED_END:
            trim_command.append("--paired")
        trim_command.extend(["--output_dir", str(trim_tmp), *raw_reads])
        trim_artifacts = []
        for canonical_read in (
            resolved.trimmed_fastq.r1,
            resolved.trimmed_fastq.r2,
        ):
            if canonical_read is None:
                continue
            trim_artifacts.append(
                artifact_contract(
                    "fastq",
                    trim_tmp / "canonical" / canonical_read.name,
                    canonical_read,
                    validator="gzip_fastq",
                )
            )
        trimming_report = output / "qc" / "trimming" / sample_id / "trimming_report.txt"
        trim_artifacts.append(
            artifact_contract(
                "text",
                trim_tmp / "canonical" / "trimming_report.txt",
                trimming_report,
                validator="text_nonempty",
            )
        )
        steps.append(
            _record(
                "trim",
                "sample",
                sample_id,
                [
                    command_action(trim_command),
                    internal_action(
                        "normalize_trim_galore_outputs",
                        output_dir=str(trim_tmp),
                        read_layout=config.read_layout.value,
                        canonical_temporary_paths=[
                            artifact["temporary_path"] for artifact in trim_artifacts
                        ],
                        write_mode="atomic_replace",
                    ),
                ],
                trim_artifacts,
            )
        )

        trimmed_qc_tmp = _step_tmp(
            output, "samples", sample_id, "fastqc_trimmed"
        )
        trimmed_qc = output / "qc" / "fastqc" / "trimmed" / sample_id
        steps.append(
            _record(
                "fastqc_trimmed",
                "sample",
                sample_id,
                [
                    command_action(
                        [
                            _tool(tools, "fastqc"),
                            "--threads",
                            str(config.threads),
                            "--outdir",
                            str(trimmed_qc_tmp),
                            *trimmed_reads,
                        ]
                    )
                ],
                _fastqc_artifacts(trimmed_reads, trimmed_qc_tmp, trimmed_qc),
            )
        )

        align_tmp = _step_tmp(output, "samples", sample_id, "align")
        temporary_collated = align_tmp / f"{sample_id}.name-collated.bam"
        canonical_collated = (
            output
            / "intermediate"
            / "alignment"
            / f"{sample_id}.name-collated.bam"
        )
        align = [_tool(tools, "bowtie2"), "--very-sensitive"]
        if config.read_layout is ReadLayout.PAIRED_END:
            align.extend(["--no-mixed", "--no-discordant"])
            if config.allow_dovetail:
                align.append("--dovetail")
            align.extend(
                [
                    "-x",
                    str(config.bowtie2_index),
                    "-1",
                    trimmed_reads[0],
                    "-2",
                    trimmed_reads[1],
                ]
            )
        else:
            align.extend(["-x", str(config.bowtie2_index), "-U", trimmed_reads[0]])
        align.extend(["-p", str(thread_allocation.bowtie2_threads)])
        steps.append(
            _record(
                "align",
                "sample",
                sample_id,
                [
                    pipeline_action(
                        [
                            align,
                            [
                                _tool(tools, "samtools"),
                                "collate",
                                "-@",
                                str(thread_allocation.collate_workers),
                                "-o",
                                str(temporary_collated),
                                "-",
                            ],
                        ]
                    )
                ],
                [
                    artifact_contract(
                        "bam",
                        temporary_collated,
                        canonical_collated,
                        validator="bam_nonempty",
                    )
                ],
            )
        )

        bam_tmp = _step_tmp(output, "samples", sample_id, "bam_process")
        coordinate_bam = bam_tmp / "coordinate.bam"
        final_bam_tmp = bam_tmp / f"{sample_id}.final.bam"
        final_bai_tmp = Path(str(final_bam_tmp) + ".bai")
        final_bam = output / "bam" / f"{sample_id}.final.bam"
        markdup_tmp = bam_tmp / f"{sample_id}.markdup.txt"
        markdup = [
            _tool(tools, "samtools"),
            "markdup",
            "-s",
            "-f",
            str(markdup_tmp),
        ]
        if config.duplicate_policy == "remove":
            markdup.append("-r")
        markdup.extend([str(coordinate_bam), "-"])
        filter_command = [
            _tool(tools, "samtools"),
            "view",
            "-b",
            "-q",
            str(config.min_mapq),
        ]
        if config.read_layout is ReadLayout.PAIRED_END:
            filter_command.extend(["-f", "2", "-F", "2828"])
        else:
            filter_command.extend(["-F", "2820"])
        filter_command.extend(["-o", str(final_bam_tmp), "-"])
        qc_paths = {
            "flagstat": bam_tmp / f"{sample_id}.flagstat.txt",
            "stats": bam_tmp / f"{sample_id}.stats.txt",
            "idxstats": bam_tmp / f"{sample_id}.idxstats.txt",
        }
        bam_actions = [
            pipeline_action(
                [
                    [
                        _tool(tools, "samtools"),
                        "fixmate",
                        "-m",
                        str(canonical_collated),
                        "-",
                    ],
                    [
                        _tool(tools, "samtools"),
                        "sort",
                        "-@",
                        str(thread_allocation.sort_workers),
                        "-o",
                        str(coordinate_bam),
                        "-",
                    ],
                ]
            ),
            pipeline_action([markdup, filter_command]),
            command_action(
                [
                    _tool(tools, "samtools"),
                    "index",
                    "-@",
                    str(thread_allocation.samtools_workers),
                    str(final_bam_tmp),
                ]
            ),
            command_action(
                [
                    _tool(tools, "samtools"),
                    "flagstat",
                    "-@",
                    str(thread_allocation.samtools_workers),
                    str(final_bam_tmp),
                ],
                stdout_path=qc_paths["flagstat"],
            ),
            command_action(
                [
                    _tool(tools, "samtools"),
                    "stats",
                    "-@",
                    str(thread_allocation.samtools_workers),
                    str(final_bam_tmp),
                ],
                stdout_path=qc_paths["stats"],
            ),
            command_action(
                [_tool(tools, "samtools"), "idxstats", str(final_bam_tmp)],
                stdout_path=qc_paths["idxstats"],
            ),
        ]
        bam_artifacts = [
            artifact_contract(
                "bam", final_bam_tmp, final_bam, validator="final_bam"
            ),
            artifact_contract(
                "bai",
                final_bai_tmp,
                Path(str(final_bam) + ".bai"),
                validator="bam_index",
            ),
            artifact_contract(
                "text",
                markdup_tmp,
                output / "qc" / "samtools" / f"{sample_id}.markdup.txt",
                validator="text_nonempty",
            ),
        ]
        for kind, temporary in qc_paths.items():
            bam_artifacts.append(
                artifact_contract(
                    "text",
                    temporary,
                    output / "qc" / "samtools" / f"{sample_id}.{kind}.txt",
                    validator="text_nonempty",
                )
            )
        steps.append(
            _record(
                "bam_process",
                "sample",
                sample_id,
                bam_actions,
                bam_artifacts,
            )
        )

        coverage_tmp = _step_tmp(output, "samples", sample_id, "coverage")
        bigwig_tmp = coverage_tmp / f"{sample_id}.normalized.bw"
        bigwig = output / "bigwig" / f"{sample_id}.normalized.bw"
        coverage = [
            _tool(tools, "bamCoverage"),
            "--bam",
            str(final_bam),
            "--outFileName",
            str(bigwig_tmp),
            "--binSize",
            str(config.bin_size),
            "--numberOfProcessors",
            str(config.threads),
            "--normalizeUsing",
            config.normalization,
        ]
        if config.normalization == "RPGC":
            coverage.extend(
                ["--effectiveGenomeSize", str(config.effective_genome_size)]
            )
        steps.append(
            _record(
                "coverage",
                "sample",
                sample_id,
                [command_action(coverage)],
                [
                    artifact_contract(
                        "bigwig", bigwig_tmp, bigwig, validator="bigwig_magic"
                    )
                ],
            )
        )

        if includes("peak") and sample.peak_type != "input":
            extension = (
                "broadPeak" if sample.peak_type == "broad" else "narrowPeak"
            )
            peak_root = output / "peaks" / "replicate" / sample_id
            peak_name = f"{sample_id}.peaks.{extension}"
            control_bam = (
                None
                if sample.control == "."
                else output / "bam" / f"{sample.control}.final.bam"
            )
            actions, artifacts = _peak_plan(
                config=config,
                tools=tools,
                scope="sample",
                scope_id=sample_id,
                name=sample_id,
                treatment_bam=final_bam,
                control_bam=control_bam,
                peak_type=sample.peak_type,
                qvalue=sample.qvalue,
                canonical_root=peak_root,
                canonical_name=peak_name,
            )
            steps.append(
                _record("peak", "sample", sample_id, actions, artifacts)
            )

    if includes("multiqc"):
        multiqc_tmp = _step_tmp(output, "pipeline", "run", "multiqc") / "output"
        steps.append(
            _record(
                "multiqc",
                "pipeline",
                "run",
                [
                    command_action(
                        [
                            _tool(tools, "multiqc"),
                            str(output),
                            "--outdir",
                            str(multiqc_tmp),
                        ]
                    )
                ],
                [
                    artifact_contract(
                        "directory",
                        multiqc_tmp,
                        output / "multiqc",
                        validator="multiqc_output",
                    )
                ],
            )
        )

    for group in groups if includes("pool_bam") else ():
        pool_tmp = _step_tmp(output, "groups", group.group_id, "pool_bam")
        pooled_bam_tmp = pool_tmp / f"{group.group_id}.pooled.bam"
        pooled_bai_tmp = Path(str(pooled_bam_tmp) + ".bai")
        pooled_bam = output / "bam" / "pooled" / f"{group.group_id}.pooled.bam"
        member_bams = [
            str(output / "bam" / f"{sample_id}.final.bam")
            for sample_id in group.sample_ids
        ]
        steps.append(
            _record(
                "pool_bam",
                "group",
                group.group_id,
                [
                    command_action(
                        [
                            _tool(tools, "samtools"),
                            "merge",
                            "-@",
                            str(thread_allocation.samtools_workers),
                            str(pooled_bam_tmp),
                            *member_bams,
                        ]
                    ),
                    command_action(
                        [
                            _tool(tools, "samtools"),
                            "index",
                            "-@",
                            str(thread_allocation.samtools_workers),
                            str(pooled_bam_tmp),
                        ]
                    ),
                ],
                [
                    artifact_contract(
                        "bam", pooled_bam_tmp, pooled_bam, validator="pooled_bam"
                    ),
                    artifact_contract(
                        "bai",
                        pooled_bai_tmp,
                        Path(str(pooled_bam) + ".bai"),
                        validator="bam_index",
                    ),
                ],
            )
        )

        if not includes("pooled_peak"):
            continue

        extension = "broadPeak" if group.peak_type == "broad" else "narrowPeak"
        pooled_root = output / "peaks" / "pooled" / group.group_id
        pooled_name = f"{group.group_id}.pooled.peaks.{extension}"
        group_control = (
            None
            if group.control == "."
            else output / "bam" / f"{group.control}.final.bam"
        )
        pooled_actions, pooled_artifacts = _peak_plan(
            config=config,
            tools=tools,
            scope="group",
            scope_id=group.group_id,
            name=group.group_id,
            treatment_bam=pooled_bam,
            control_bam=group_control,
            peak_type=group.peak_type,
            qvalue=group.qvalue,
            canonical_root=pooled_root,
            canonical_name=pooled_name,
        )
        steps.append(
            _record(
                "pooled_peak",
                "group",
                group.group_id,
                pooled_actions,
                pooled_artifacts,
            )
        )

        if not includes("consensus_peak"):
            continue

        replicate_peaks = tuple(
            output
            / "peaks"
            / "replicate"
            / sample_id
            / f"{sample_id}.peaks.{extension}"
            for sample_id in group.sample_ids
        )
        consensus_actions, consensus_artifacts, consensus_metadata = (
            build_consensus_plan(
                bedtools=_tool(tools, "bedtools"),
                output_dir=output,
                group_id=group.group_id,
                sample_ids=group.sample_ids,
                replicate_peaks=replicate_peaks,
            )
        )
        steps.append(
            _record(
                "consensus_peak",
                "group",
                group.group_id,
                consensus_actions,
                consensus_artifacts,
                metadata=consensus_metadata,
            )
        )

    if not includes("report"):
        steps = [
            step
            for step in steps
            if PHASES.index(str(step["phase"])) <= stop_index
        ]
        return {
            "schema_version": 2,
            "dry_run": config.dry_run,
            "read_layout": config.read_layout.value,
            "stop_after": selected_stop,
            "thread_allocation": thread_allocation.as_plan_metadata(),
            "artifact_contract": {
                "command_outputs": "write temporary_path only",
                "validation": "run the declared validator on temporary_path",
                "promotion": (
                    "validated files use atomic replacement; validated directories "
                    "use a same-parent staged replacement"
                ),
                "failure": "never promote partial or invalid output",
            },
            "steps": _phase_major(steps),
        }

    report_tmp = _step_tmp(output, "pipeline", "run", "report")
    report_files = {
        "run_summary.tsv": "summary_tsv",
        "run_summary.json": "summary_json",
        "warnings.tsv": "warnings_tsv",
    }
    report_actions = [
        internal_action(
            "build_run_summary",
            temporary_output_dir=str(report_tmp),
            output_files=list(report_files),
            output_manifest_path=str(report_tmp / "output_manifest.json"),
            write_mode="atomic_replace",
        ),
    ]
    report_artifacts = []
    for filename, validator in report_files.items():
        report_artifacts.append(
            artifact_contract(
                "report",
                report_tmp / filename,
                output / "report" / filename,
                validator=validator,
            )
        )
    report_artifacts.append(
        artifact_contract(
            "json",
            report_tmp / "output_manifest.json",
            output / "output_manifest.json",
            validator="output_manifest",
        )
    )
    steps.append(
        _record(
            "report",
            "pipeline",
            "run",
            report_actions,
            report_artifacts,
        )
    )
    return {
        "schema_version": 2,
        "dry_run": config.dry_run,
        "read_layout": config.read_layout.value,
        "stop_after": selected_stop,
        "thread_allocation": thread_allocation.as_plan_metadata(),
        "artifact_contract": {
            "command_outputs": "write temporary_path only",
            "validation": "run the declared validator on temporary_path",
            "promotion": (
                "validated files use atomic replacement; validated directories "
                "use a same-parent staged replacement"
            ),
            "failure": "never promote partial or invalid output",
        },
        "steps": _phase_major(steps),
    }
