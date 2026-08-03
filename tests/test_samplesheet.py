from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from wulfpeak.samplesheet import (
    REQUIRED_COLUMNS,
    SamplesheetValidationError,
    load_samplesheet,
)


def write_sheet(path: Path, rows: list[dict[str, str]], columns=REQUIRED_COLUMNS) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def valid_rows() -> list[dict[str, str]]:
    return [
        {
            "input_id": "chip_S1",
            "sample_id": "chip_1",
            "group_id": "chip",
            "peak_type": "narrow",
            "control": "input_1",
            "qvalue": "0.01",
        },
        {
            "input_id": "input_S2",
            "sample_id": "input_1",
            "group_id": "input",
            "peak_type": "input",
            "control": ".",
            "qvalue": ".",
        },
    ]


class SamplesheetTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.sheet = Path(self.temporary.name) / "samples.tsv"

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_loads_valid_sheet(self) -> None:
        write_sheet(self.sheet, valid_rows())
        samples = load_samplesheet(self.sheet)
        self.assertEqual([sample.sample_id for sample in samples], ["chip_1", "input_1"])

    def test_column_order_is_not_significant(self) -> None:
        write_sheet(self.sheet, valid_rows(), tuple(reversed(REQUIRED_COLUMNS)))
        self.assertEqual(len(load_samplesheet(self.sheet)), 2)

    def test_requires_exact_column_set(self) -> None:
        cases = (
            (REQUIRED_COLUMNS[:-1], "missing required columns: qvalue"),
            (REQUIRED_COLUMNS + ("fastq1",), "unsupported columns in v0.1.0: fastq1"),
        )
        for columns, message in cases:
            with self.subTest(columns=columns):
                rows = [
                    {key: value for key, value in row.items() if key in columns}
                    for row in valid_rows()
                ]
                for row in rows:
                    if "fastq1" in columns:
                        row["fastq1"] = "legacy.fastq.gz"
                write_sheet(self.sheet, rows, columns)
                with self.assertRaisesRegex(SamplesheetValidationError, message):
                    load_samplesheet(self.sheet)

    def test_rejects_duplicate_columns(self) -> None:
        duplicate_columns = REQUIRED_COLUMNS + ("sample_id",)
        self.sheet.write_text(
            "\t".join(duplicate_columns) + "\n"
            "raw\tsample\tgroup\tnarrow\t.\t0.01\tduplicate\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(
            SamplesheetValidationError, "columns must not be repeated: sample_id"
        ):
            load_samplesheet(self.sheet)

    def test_rejects_rows_with_extra_fields(self) -> None:
        self.sheet.write_text(
            "\t".join(REQUIRED_COLUMNS) + "\n"
            "raw\tsample\tgroup\tnarrow\t.\t0.01\tunexpected\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(
            SamplesheetValidationError,
            "contains more tab-separated fields than the header",
        ):
            load_samplesheet(self.sheet)

    def test_required_identifiers_must_be_nonempty(self) -> None:
        for field in ("input_id", "sample_id", "group_id"):
            with self.subTest(field=field):
                rows = valid_rows()
                rows[0][field] = ""
                write_sheet(self.sheet, rows)
                with self.assertRaisesRegex(
                    SamplesheetValidationError, f"{field} must be non-empty"
                ):
                    load_samplesheet(self.sheet)

    def test_ids_must_be_unique(self) -> None:
        for field in ("input_id", "sample_id"):
            with self.subTest(field=field):
                rows = valid_rows()
                rows[1][field] = rows[0][field]
                write_sheet(self.sheet, rows)
                with self.assertRaisesRegex(
                    SamplesheetValidationError, f"{field} values must be unique"
                ):
                    load_samplesheet(self.sheet)

    def test_peak_type_must_be_supported(self) -> None:
        rows = valid_rows()
        rows[0]["peak_type"] = "sharp"
        write_sheet(self.sheet, rows)
        with self.assertRaisesRegex(
            SamplesheetValidationError, "peak_type must be narrow, broad, or input"
        ):
            load_samplesheet(self.sheet)

    def test_peak_qvalue_must_be_finite_numeric(self) -> None:
        for bad_qvalue in (".", "abc", "NaN", "inf"):
            with self.subTest(bad_qvalue=bad_qvalue):
                rows = valid_rows()
                rows[0]["qvalue"] = bad_qvalue
                write_sheet(self.sheet, rows)
                with self.assertRaisesRegex(
                    SamplesheetValidationError, "qvalue must be numeric and finite"
                ):
                    load_samplesheet(self.sheet)

    def test_peak_qvalue_must_be_in_macs_range(self) -> None:
        for bad_qvalue in ("-1", "0", "1.1"):
            with self.subTest(bad_qvalue=bad_qvalue):
                rows = valid_rows()
                rows[0]["qvalue"] = bad_qvalue
                write_sheet(self.sheet, rows)
                with self.assertRaisesRegex(
                    SamplesheetValidationError,
                    "qvalue must satisfy 0 < qvalue <= 1",
                ):
                    load_samplesheet(self.sheet)

    def test_peak_qvalue_accepts_range_boundaries(self) -> None:
        for valid_qvalue in ("1", "0.01"):
            with self.subTest(valid_qvalue=valid_qvalue):
                rows = valid_rows()
                rows[0]["qvalue"] = valid_qvalue
                write_sheet(self.sheet, rows)
                self.assertEqual(load_samplesheet(self.sheet)[0].qvalue, valid_qvalue)

    def test_non_input_group_requires_one_peak_type(self) -> None:
        rows = valid_rows()
        rows.insert(
            1,
            {
                "input_id": "chip_S2",
                "sample_id": "chip_2",
                "group_id": "chip",
                "peak_type": "broad",
                "control": "input_1",
                "qvalue": "0.01",
            },
        )
        write_sheet(self.sheet, rows)

        with self.assertRaisesRegex(
            SamplesheetValidationError,
            "group_id 'chip': non-input samples must share the same peak_type; "
            "found broad, narrow",
        ):
            load_samplesheet(self.sheet)

    def test_input_samples_do_not_participate_in_group_peak_type_validation(self) -> None:
        rows = valid_rows()
        rows[1]["group_id"] = "chip"
        write_sheet(self.sheet, rows)

        self.assertEqual(len(load_samplesheet(self.sheet)), 2)

    def test_non_input_group_requires_one_qvalue(self) -> None:
        rows = valid_rows()
        rows.insert(
            1,
            {
                "input_id": "chip_S2",
                "sample_id": "chip_2",
                "group_id": "chip",
                "peak_type": "narrow",
                "control": "input_1",
                "qvalue": "0.05",
            },
        )
        write_sheet(self.sheet, rows)

        with self.assertRaisesRegex(
            SamplesheetValidationError,
            "group_id 'chip': non-input samples must use the same qvalue; "
            "found 0.01, 0.05",
        ):
            load_samplesheet(self.sheet)

    def test_group_qvalue_comparison_is_numeric(self) -> None:
        rows = valid_rows()
        rows.insert(
            1,
            {
                "input_id": "chip_S2",
                "sample_id": "chip_2",
                "group_id": "chip",
                "peak_type": "narrow",
                "control": "input_1",
                "qvalue": "1e-2",
            },
        )
        write_sheet(self.sheet, rows)

        self.assertEqual(len(load_samplesheet(self.sheet)), 3)

    def test_input_row_requires_dots(self) -> None:
        for field, value in (("qvalue", "0.01"), ("control", "chip_1")):
            with self.subTest(field=field):
                rows = valid_rows()
                rows[1][field] = value
                write_sheet(self.sheet, rows)
                with self.assertRaisesRegex(
                    SamplesheetValidationError, f"input samples must use {field}='\\.'"
                ):
                    load_samplesheet(self.sheet)

    def test_peak_control_must_reference_input_sample(self) -> None:
        rows = valid_rows()
        rows[0]["control"] = "unknown"
        write_sheet(self.sheet, rows)
        with self.assertRaisesRegex(SamplesheetValidationError, "must match an existing sample_id"):
            load_samplesheet(self.sheet)

    def test_peak_control_may_be_dot(self) -> None:
        rows = valid_rows()
        rows[0]["control"] = "."
        write_sheet(self.sheet, rows)
        self.assertEqual(len(load_samplesheet(self.sheet)), 2)


if __name__ == "__main__":
    unittest.main()
