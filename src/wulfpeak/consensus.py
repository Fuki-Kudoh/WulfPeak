"""Strict-majority consensus planning without executing analysis commands."""

from __future__ import annotations

from pathlib import Path

from .plan_schema import artifact_contract, command_action, internal_action


SUPPORT_COLUMNS = (
    "chrom",
    "start",
    "end",
    "support_count",
    "required_count",
    "replicate_count",
    "supporting_samples",
)


def strict_majority_threshold(replicate_count: int) -> int:
    if replicate_count <= 0:
        raise ValueError("replicate_count must be positive")
    return replicate_count // 2 + 1


def build_consensus_plan(
    *,
    bedtools: str,
    output_dir: Path,
    group_id: str,
    sample_ids: tuple[str, ...],
    replicate_peaks: tuple[Path, ...],
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    if len(sample_ids) != len(replicate_peaks) or not sample_ids:
        raise ValueError("consensus inputs must contain one peak per sample")
    step_tmp = (
        output_dir
        / "intermediate"
        / ".steps"
        / "groups"
        / group_id
        / "consensus_peak"
    )
    actions: list[dict[str, object]] = []
    normalized: list[Path] = []
    normalization: list[dict[str, str]] = []
    for sample_id, peak in zip(sample_ids, replicate_peaks):
        sorted_peak = step_tmp / "sorted" / f"{sample_id}.bed"
        normalized_peak = step_tmp / "normalized" / f"{sample_id}.bed"
        actions.append(
            command_action(
                [bedtools, "sort", "-i", str(peak)], stdout_path=sorted_peak
            )
        )
        actions.append(
            command_action(
                [bedtools, "merge", "-i", str(sorted_peak)],
                stdout_path=normalized_peak,
            )
        )
        normalized.append(normalized_peak)
        normalization.append(
            {
                "sample_id": sample_id,
                "source_peak": str(peak),
                "normalized_bed3": str(normalized_peak),
            }
        )

    multiinter = step_tmp / "multiinter.tsv"
    actions.append(
        command_action(
            [
                bedtools,
                "multiinter",
                "-header",
                "-names",
                *sample_ids,
                "-i",
                *(str(path) for path in normalized),
            ],
            stdout_path=multiinter,
        )
    )
    required = strict_majority_threshold(len(sample_ids))
    qualifying = step_tmp / "qualifying.bed"
    support_tmp = step_tmp / f"{group_id}.consensus.support.tsv"
    actions.append(
        internal_action(
            "filter_consensus_support",
            input_path=str(multiinter),
            minimum_support=required,
            replicate_count=len(sample_ids),
            sample_ids=list(sample_ids),
            supporting_samples_source="bedtools_multiinter_names",
            qualifying_bed3_path=str(qualifying),
            support_tsv_path=str(support_tmp),
            support_columns=list(SUPPORT_COLUMNS),
            write_mode="atomic_replace",
        )
    )
    consensus_tmp = step_tmp / f"{group_id}.consensus.bed"
    actions.append(
        command_action(
            [bedtools, "merge", "-i", str(qualifying)],
            stdout_path=consensus_tmp,
        )
    )
    consensus = (
        output_dir
        / "peaks"
        / "consensus"
        / group_id
        / f"{group_id}.consensus.bed"
    )
    support = consensus.with_name(f"{group_id}.consensus.support.tsv")
    artifacts = [
        artifact_contract(
            "bed3", consensus_tmp, consensus, validator="bed3_sorted_or_empty"
        ),
        artifact_contract(
            "tsv", support_tmp, support, validator="consensus_support_tsv"
        ),
    ]
    metadata = {
        "algorithm": "strict_majority",
        "replicate_count": len(sample_ids),
        "required_count": required,
        "sample_ids": list(sample_ids),
        "normalization": normalization,
        "support_columns": list(SUPPORT_COLUMNS),
    }
    return actions, artifacts, metadata
