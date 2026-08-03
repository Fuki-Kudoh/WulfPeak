"""Read and validate the compact WulfPeak v0.1.0 samplesheet."""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path


REQUIRED_COLUMNS = (
    "input_id",
    "sample_id",
    "group_id",
    "peak_type",
    "control",
    "qvalue",
)
PEAK_TYPES = frozenset({"narrow", "broad", "input"})


class SamplesheetValidationError(ValueError):
    """Raised when a samplesheet violates the v0.1.0 contract."""

    def __init__(self, errors: list[str]):
        self.errors = errors
        super().__init__("Samplesheet validation failed:\n" + "\n".join(
            f"  - {error}" for error in errors
        ))


@dataclass(frozen=True)
class Sample:
    input_id: str
    sample_id: str
    group_id: str
    peak_type: str
    control: str
    qvalue: str


def _duplicates(values: list[str]) -> list[str]:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    return sorted(duplicates)


def _is_numeric(value: str) -> bool:
    try:
        return math.isfinite(float(value))
    except ValueError:
        return False


def load_samplesheet(path: str | Path) -> list[Sample]:
    """Load a tab-separated samplesheet and enforce the v0.1.0 schema."""

    samplesheet = Path(path)
    errors: list[str] = []
    try:
        handle = samplesheet.open(newline="", encoding="utf-8")
    except OSError as exc:
        raise SamplesheetValidationError([f"cannot read {samplesheet}: {exc}"]) from exc

    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        actual_columns = tuple(reader.fieldnames or ())
        missing = [name for name in REQUIRED_COLUMNS if name not in actual_columns]
        extra = [name for name in actual_columns if name not in REQUIRED_COLUMNS]
        duplicate_columns = _duplicates(list(actual_columns))
        if missing:
            errors.append(f"missing required columns: {', '.join(missing)}")
        if extra:
            errors.append(
                "unsupported columns in v0.1.0: " + ", ".join(extra)
            )
        if duplicate_columns:
            errors.append(
                "columns must not be repeated: " + ", ".join(duplicate_columns)
            )
        if errors:
            raise SamplesheetValidationError(errors)

        rows = list(reader)

    if not rows:
        raise SamplesheetValidationError(["samplesheet contains no sample rows"])

    samples: list[Sample] = []
    for row_number, row in enumerate(rows, start=2):
        if None in row:
            errors.append(
                f"row {row_number}: contains more tab-separated fields than the header"
            )
        values = {name: (row.get(name) or "").strip() for name in REQUIRED_COLUMNS}
        for name in REQUIRED_COLUMNS:
            if not values[name]:
                errors.append(f"row {row_number}: {name} must be non-empty")

        peak_type = values["peak_type"]
        control = values["control"]
        qvalue = values["qvalue"]
        if peak_type and peak_type not in PEAK_TYPES:
            errors.append(
                f"row {row_number}: peak_type must be narrow, broad, or input; "
                f"got {peak_type!r}"
            )
        elif peak_type == "input":
            if control != ".":
                errors.append(f"row {row_number}: input samples must use control='.'")
            if qvalue != ".":
                errors.append(f"row {row_number}: input samples must use qvalue='.'")
        elif peak_type in {"narrow", "broad"}:
            if not _is_numeric(qvalue):
                errors.append(
                    f"row {row_number}: qvalue must be numeric for {peak_type} samples; "
                    f"got {qvalue!r}"
                )

        samples.append(Sample(**values))

    duplicate_sample_ids = _duplicates([sample.sample_id for sample in samples])
    if duplicate_sample_ids:
        errors.append(
            "sample_id values must be unique; duplicates: "
            + ", ".join(duplicate_sample_ids)
        )

    duplicate_input_ids = _duplicates([sample.input_id for sample in samples])
    if duplicate_input_ids:
        errors.append(
            "input_id values must be unique; duplicates: "
            + ", ".join(duplicate_input_ids)
        )

    input_sample_ids = {
        sample.sample_id for sample in samples if sample.peak_type == "input"
    }
    for row_number, sample in enumerate(samples, start=2):
        if sample.peak_type in {"narrow", "broad"} and sample.control != ".":
            if sample.control not in input_sample_ids:
                errors.append(
                    f"row {row_number}: control {sample.control!r} must match an "
                    "existing sample_id with peak_type=input, or be '.'"
                )

    if errors:
        raise SamplesheetValidationError(errors)
    return samples
