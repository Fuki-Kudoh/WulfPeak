"""Centralized per-sample resource allocation policies."""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class ThreadAllocation:
    """Thread arguments derived from one approximate per-sample CPU budget."""

    budget: int
    bowtie2_threads: int
    collate_workers: int
    sort_workers: int
    samtools_workers: int

    def as_plan_metadata(self) -> dict[str, int]:
        return asdict(self)


def allocate_threads(budget: int) -> ThreadAllocation:
    """Split a sample budget without giving it independently to piped processes.

    Samtools ``-@`` counts additional workers, not the command's main thread.
    A streaming pipeline necessarily has two main processes, so budget 1 uses
    no additional workers but may briefly consume two CPUs. Budgets of 2 or
    more stay within the requested total for the concurrent pipelines.
    """

    if budget <= 0:
        raise ValueError("thread budget must be positive")

    collate_total = max(1, budget // 4)
    bowtie2_threads = max(1, budget - collate_total)
    return ThreadAllocation(
        budget=budget,
        bowtie2_threads=bowtie2_threads,
        collate_workers=collate_total - 1,
        sort_workers=max(0, budget - 2),
        samtools_workers=max(0, budget - 1),
    )
