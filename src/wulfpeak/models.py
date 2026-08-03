"""Core immutable models shared across CLI and pipeline layers."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path


class ReadLayout(str, Enum):
    PAIRED_END = "paired_end"
    SINGLE_END = "single_end"


class Assay(str, Enum):
    CHIPSEQ = "chipseq"
    CUTANDTAG = "cutandtag"


@dataclass(frozen=True)
class FastqInput:
    r1: Path
    r2: Path | None

    @property
    def layout(self) -> ReadLayout:
        if self.r2 is None:
            return ReadLayout.SINGLE_END
        return ReadLayout.PAIRED_END


@dataclass(frozen=True)
class ReplicateGroup:
    group_id: str
    sample_ids: tuple[str, ...]
    peak_type: str
    control: str
    qvalue: str


@dataclass(frozen=True)
class ToolInfo:
    name: str
    path: Path
    version: str | None
    version_command: tuple[str, ...]


@dataclass(frozen=True)
class ValidationResult:
    valid: bool
    checks: dict[str, object]
    warnings: tuple[str, ...] = ()
