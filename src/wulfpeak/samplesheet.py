"""Read and validate the compact WulfPeak v0.1.0 samplesheet."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
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
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


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


def _parse_qvalue(value: str) -> Decimal | None:
    try:
        qvalue = Decimal(value)
    except InvalidOperation:
        return None
    return qvalue if qvalue.is_finite() else None


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
        raw_values = {name: (row.get(name) or "") for name in REQUIRED_COLUMNS}
        values = {name: value.strip() for name, value in raw_values.items()}
        for name in REQUIRED_COLUMNS:
            if not values[name]:
                errors.append(f"row {row_number}: {name} must be non-empty")

        for name in ("sample_id", "group_id"):
            value = values[name]
            if raw_values[name] != value:
                errors.append(
                    f"row {row_number}: {name} must not have leading or trailing whitespace"
                )
            if value in {".", ".."} or (value and not SAFE_ID.fullmatch(value)):
                errors.append(
                    f"row {row_number}: {name} must match {SAFE_ID.pattern!r} and "
                    "must not be '.' or '..'"
                )

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
            parsed_qvalue = _parse_qvalue(qvalue)
            if parsed_qvalue is None:
                errors.append(
                    f"row {row_number}: qvalue must be numeric and finite for "
                    f"{peak_type} samples; "
                    f"got {qvalue!r}"
                )
            elif not Decimal("0") < parsed_qvalue <= Decimal("1"):
                errors.append(
                    f"row {row_number}: qvalue must satisfy 0 < qvalue <= 1 for "
                    f"{peak_type} samples; got {qvalue!r}"
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

    non_input_groups: dict[str, list[Sample]] = {}
    for sample in samples:
        if sample.peak_type in {"narrow", "broad"}:
            non_input_groups.setdefault(sample.group_id, []).append(sample)

    for group_id, group_samples in non_input_groups.items():
        peak_types = sorted({sample.peak_type for sample in group_samples})
        if len(peak_types) > 1:
            errors.append(
                f"group_id {group_id!r}: non-input samples must share the same "
                f"peak_type; found {', '.join(peak_types)}"
            )

        qvalues: dict[Decimal, set[str]] = {}
        for sample in group_samples:
            parsed_qvalue = _parse_qvalue(sample.qvalue)
            if (
                parsed_qvalue is not None
                and Decimal("0") < parsed_qvalue <= Decimal("1")
            ):
                qvalues.setdefault(parsed_qvalue, set()).add(sample.qvalue)
        if len(qvalues) > 1:
            displayed_qvalues = sorted(
                value
                for equivalent_values in qvalues.values()
                for value in equivalent_values
            )
            errors.append(
                f"group_id {group_id!r}: non-input samples must use the same "
                f"qvalue; found {', '.join(displayed_qvalues)}"
            )

        controls = sorted({sample.control for sample in group_samples})
        if len(controls) > 1:
            errors.append(
                f"group_id {group_id!r}: non-input samples must use the same "
                f"control; found {', '.join(controls)}"
            )

    if errors:
        raise SamplesheetValidationError(errors)
    return samples
