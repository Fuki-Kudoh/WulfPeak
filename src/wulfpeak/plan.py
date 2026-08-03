"""Deterministic command-plan construction for the execution-foundation PR."""

from __future__ import annotations

from pathlib import Path

from .config import RunConfig
from .manifest import ResolvedSample
from .models import ReadLayout, ReplicateGroup, ToolInfo


def _tool(tools: dict[str, ToolInfo], name: str) -> str:
    return str(tools[name].path)


def _step_tmp(output: Path, scope: str, scope_id: str, phase: str) -> str:
    return str(output / "intermediate" / ".steps" / scope / scope_id / phase)


def _record(
    phase: str,
    scope: str,
    scope_id: str,
    commands: list[list[str]],
    outputs: list[str],
    *,
    pipeline: bool = False,
) -> dict[str, object]:
    return {
        "phase": phase,
        "scope": scope,
        "scope_id": scope_id,
        "commands": commands,
        "outputs": outputs,
        "execution": "pipeline" if pipeline else "sequential",
    }


def build_command_plan(
    config: RunConfig,
    samples: tuple[ResolvedSample, ...],
    groups: tuple[ReplicateGroup, ...],
    tools: dict[str, ToolInfo],
) -> dict[str, object]:
    """Return analysis commands without executing any of them."""

    output = config.output_dir
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

        raw_qc_dir = output / "qc" / "fastqc" / "raw" / sample_id
        steps.append(
            _record(
                "fastqc_raw",
                "sample",
                sample_id,
                [
                    [
                        _tool(tools, "fastqc"),
                        "--threads",
                        str(config.threads),
                        "--outdir",
                        _step_tmp(output, "samples", sample_id, "fastqc_raw"),
                        *raw_reads,
                    ]
                ],
                [str(raw_qc_dir)],
            )
        )
        trim_command = [
            _tool(tools, "trim_galore"),
            "--cores",
            str(config.threads),
        ]
        if config.read_layout is ReadLayout.PAIRED_END:
            trim_command.append("--paired")
        trim_command.extend(
            ["--output_dir", _step_tmp(output, "samples", sample_id, "trim"), *raw_reads]
        )
        steps.append(
            _record(
                "trim",
                "sample",
                sample_id,
                [trim_command],
                trimmed_reads,
            )
        )
        steps.append(
            _record(
                "fastqc_trimmed",
                "sample",
                sample_id,
                [
                    [
                        _tool(tools, "fastqc"),
                        "--threads",
                        str(config.threads),
                        "--outdir",
                        _step_tmp(output, "samples", sample_id, "fastqc_trimmed"),
                        *trimmed_reads,
                    ]
                ],
                [str(output / "qc" / "fastqc" / "trimmed" / sample_id)],
            )
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
        align.extend(["-p", str(config.threads)])
        unsorted_bam = (
            output / "intermediate" / "alignment" / f"{sample_id}.unsorted.bam"
        )
        steps.append(
            _record(
                "align",
                "sample",
                sample_id,
                [
                    align,
                    [
                        _tool(tools, "samtools"),
                        "view",
                        "-b",
                        "-o",
                        str(unsorted_bam),
                        "-",
                    ],
                ],
                [str(unsorted_bam)],
                pipeline=True,
            )
        )

        final_bam = output / "bam" / f"{sample_id}.final.bam"
        bam_tmp = Path(_step_tmp(output, "samples", sample_id, "bam_process"))
        collated_bam = bam_tmp / "name-collated.bam"
        fixmate_bam = bam_tmp / "fixmate.bam"
        coordinate_bam = bam_tmp / "coordinate.bam"
        marked_bam = bam_tmp / "marked.bam"
        markdup_stats = output / "qc" / "samtools" / f"{sample_id}.markdup.txt"
        markdup = [
            _tool(tools, "samtools"),
            "markdup",
            "-s",
            "-f",
            str(markdup_stats),
        ]
        if config.duplicate_policy == "remove":
            markdup.append("-r")
        markdup.extend([str(coordinate_bam), str(marked_bam)])
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
        filter_command.extend(["-o", str(final_bam), str(marked_bam)])
        steps.append(
            _record(
                "bam_process",
                "sample",
                sample_id,
                [
                    [
                        _tool(tools, "samtools"),
                        "collate",
                        "-@",
                        str(config.threads),
                        "-o",
                        str(collated_bam),
                        str(unsorted_bam),
                    ],
                    [
                        _tool(tools, "samtools"),
                        "fixmate",
                        "-m",
                        str(collated_bam),
                        str(fixmate_bam),
                    ],
                    [
                        _tool(tools, "samtools"),
                        "sort",
                        "-@",
                        str(config.threads),
                        "-o",
                        str(coordinate_bam),
                        str(fixmate_bam),
                    ],
                    markdup,
                    filter_command,
                    [_tool(tools, "samtools"), "index", "-@", str(config.threads), str(final_bam)],
                    [
                        _tool(tools, "samtools"),
                        "flagstat",
                        "-@",
                        str(config.threads),
                        str(final_bam),
                    ],
                    [
                        _tool(tools, "samtools"),
                        "stats",
                        "-@",
                        str(config.threads),
                        str(final_bam),
                    ],
                    [_tool(tools, "samtools"), "idxstats", str(final_bam)],
                ],
                [
                    str(final_bam),
                    str(final_bam) + ".bai",
                    str(output / "qc" / "samtools" / f"{sample_id}.flagstat.txt"),
                    str(output / "qc" / "samtools" / f"{sample_id}.stats.txt"),
                    str(output / "qc" / "samtools" / f"{sample_id}.idxstats.txt"),
                    str(markdup_stats),
                ],
            )
        )

        bigwig = output / "bigwig" / f"{sample_id}.normalized.bw"
        coverage = [
            _tool(tools, "bamCoverage"),
            "--bam",
            str(final_bam),
            "--outFileName",
            str(bigwig),
            "--binSize",
            str(config.bin_size),
            "--numberOfProcessors",
            str(config.threads),
            "--normalizeUsing",
            config.normalization,
        ]
        if config.normalization == "RPGC":
            coverage.extend(["--effectiveGenomeSize", str(config.effective_genome_size)])
        steps.append(_record("coverage", "sample", sample_id, [coverage], [str(bigwig)]))

        if sample.peak_type != "input":
            extension = "broadPeak" if sample.peak_type == "broad" else "narrowPeak"
            peak_output = (
                output
                / "peaks"
                / "replicate"
                / sample_id
                / f"{sample_id}.peaks.{extension}"
            )
            peak = [
                _tool(tools, "macs3"),
                "callpeak",
                "-t",
                str(final_bam),
            ]
            if sample.control != ".":
                control_bam = output / "bam" / f"{sample.control}.final.bam"
                peak.extend(["-c", str(control_bam)])
            peak.extend(
                [
                    "-f",
                    "BAMPE" if config.read_layout is ReadLayout.PAIRED_END else "BAM",
                    "-g",
                    str(config.effective_genome_size),
                    "-n",
                    sample_id,
                    "-q",
                    sample.qvalue,
                    "--outdir",
                    _step_tmp(output, "samples", sample_id, "peak"),
                ]
            )
            if sample.peak_type == "broad":
                peak.extend(["--broad", "--broad-cutoff", str(config.broad_cutoff)])
            commands = [peak]
            if config.blacklist is not None:
                raw_peak = (
                    Path(_step_tmp(output, "samples", sample_id, "peak"))
                    / f"{sample_id}_peaks.{extension}"
                )
                commands.append(
                    [
                        _tool(tools, "bedtools"),
                        "intersect",
                        "-v",
                        "-a",
                        str(raw_peak),
                        "-b",
                        str(config.blacklist),
                    ]
                )
            steps.append(_record("peak", "sample", sample_id, commands, [str(peak_output)]))

    for group in groups:
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
                    [
                        _tool(tools, "samtools"),
                        "merge",
                        "-@",
                        str(config.threads),
                        str(pooled_bam),
                        *member_bams,
                    ],
                    [
                        _tool(tools, "samtools"),
                        "index",
                        "-@",
                        str(config.threads),
                        str(pooled_bam),
                    ],
                ],
                [str(pooled_bam), str(pooled_bam) + ".bai"],
            )
        )
        extension = "broadPeak" if group.peak_type == "broad" else "narrowPeak"
        pooled_peak_path = (
            output
            / "peaks"
            / "pooled"
            / group.group_id
            / f"{group.group_id}.pooled.peaks.{extension}"
        )
        pooled_peak = [_tool(tools, "macs3"), "callpeak", "-t", str(pooled_bam)]
        if group.control != ".":
            pooled_peak.extend(["-c", str(output / "bam" / f"{group.control}.final.bam")])
        pooled_peak.extend(
            [
                "-f",
                "BAMPE" if config.read_layout is ReadLayout.PAIRED_END else "BAM",
                "-g",
                str(config.effective_genome_size),
                "-n",
                group.group_id,
                "-q",
                group.qvalue,
                "--outdir",
                _step_tmp(output, "groups", group.group_id, "pooled_peak"),
            ]
        )
        if group.peak_type == "broad":
            pooled_peak.extend(["--broad", "--broad-cutoff", str(config.broad_cutoff)])
        pooled_commands = [pooled_peak]
        if config.blacklist is not None:
            raw_pooled_peak = (
                Path(_step_tmp(output, "groups", group.group_id, "pooled_peak"))
                / f"{group.group_id}_peaks.{extension}"
            )
            pooled_commands.append(
                [
                    _tool(tools, "bedtools"),
                    "intersect",
                    "-v",
                    "-a",
                    str(raw_pooled_peak),
                    "-b",
                    str(config.blacklist),
                ]
            )
        steps.append(
            _record(
                "pooled_peak",
                "group",
                group.group_id,
                pooled_commands,
                [str(pooled_peak_path)],
            )
        )

        replicate_peaks = [
            str(output / "peaks" / "replicate" / sample_id / f"{sample_id}.peaks.{extension}")
            for sample_id in group.sample_ids
        ]
        consensus = (
            output
            / "peaks"
            / "consensus"
            / group.group_id
            / f"{group.group_id}.consensus.bed"
        )
        support = consensus.with_name(f"{group.group_id}.consensus.support.tsv")
        steps.append(
            _record(
                "consensus_peak",
                "group",
                group.group_id,
                [[_tool(tools, "bedtools"), "multiinter", "-i", *replicate_peaks]],
                [str(consensus), str(support)],
            )
        )

    steps.append(
        _record(
            "report",
            "pipeline",
            "run",
            [[_tool(tools, "multiqc"), str(output), "--outdir", str(output / "multiqc")]],
            [
                str(output / "multiqc" / "multiqc_report.html"),
                str(output / "report" / "run_summary.tsv"),
                str(output / "report" / "run_summary.json"),
                str(output / "report" / "warnings.tsv"),
                str(output / "output_manifest.json"),
            ],
        )
    )
    return {
        "schema_version": 1,
        "dry_run": config.dry_run,
        "read_layout": config.read_layout.value,
        "steps": steps,
    }
