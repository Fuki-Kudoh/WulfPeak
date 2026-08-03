from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from wulfpeak.fastq import (
    AmbiguousFastqError,
    FastqNotFoundError,
    FastqPairMismatchError,
    discover_fastq_input,
    discover_fastq_pair,
)
from wulfpeak.models import ReadLayout


SUPPORTED_PAIRS = (
    ("_R1.fastq.gz", "_R2.fastq.gz"),
    ("_R1.fq.gz", "_R2.fq.gz"),
    ("_R1_001.fastq.gz", "_R2_001.fastq.gz"),
    ("_R1_001.fq.gz", "_R2_001.fq.gz"),
    ("_1.fastq.gz", "_2.fastq.gz"),
    ("_1.fq.gz", "_2.fq.gz"),
)

SUPPORTED_SINGLE = (
    ".fastq.gz",
    ".fq.gz",
    "_R1.fastq.gz",
    "_R1.fq.gz",
    "_R1_001.fastq.gz",
    "_R1_001.fq.gz",
    "_1.fastq.gz",
    "_1.fq.gz",
)


class FastqDiscoveryTests(unittest.TestCase):
    def test_discovers_each_supported_single_end_name(self) -> None:
        for suffix in SUPPORTED_SINGLE:
            with self.subTest(suffix=suffix), tempfile.TemporaryDirectory() as temporary:
                fastq_dir = Path(temporary)
                expected = fastq_dir / f"library{suffix}"
                expected.touch()

                resolved = discover_fastq_input(
                    "library", fastq_dir, ReadLayout.SINGLE_END
                )

                self.assertEqual(resolved.r1, expected.resolve())
                self.assertIsNone(resolved.r2)
                self.assertEqual(resolved.layout, ReadLayout.SINGLE_END)

    def test_single_end_rejects_multiple_supported_candidates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fastq_dir = Path(temporary)
            (fastq_dir / "library.fastq.gz").touch()
            (fastq_dir / "library_R1.fastq.gz").touch()
            with self.assertRaisesRegex(AmbiguousFastqError, "will not guess"):
                discover_fastq_input("library", fastq_dir, ReadLayout.SINGLE_END)

    def test_single_end_rejects_r2_instead_of_ignoring_it(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fastq_dir = Path(temporary)
            (fastq_dir / "library.fastq.gz").touch()
            (fastq_dir / "library_R2.fastq.gz").touch()
            with self.assertRaisesRegex(FastqPairMismatchError, "will not ignore"):
                discover_fastq_input("library", fastq_dir, ReadLayout.SINGLE_END)

    def test_single_end_reports_missing_input(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(FastqNotFoundError, "No single-end FASTQ"):
                discover_fastq_input("missing", temporary, ReadLayout.SINGLE_END)

    def test_discovers_each_supported_pair(self) -> None:
        for r1_suffix, r2_suffix in SUPPORTED_PAIRS:
            with self.subTest(r1_suffix=r1_suffix, r2_suffix=r2_suffix):
                with tempfile.TemporaryDirectory() as temporary:
                    fastq_dir = Path(temporary)
                    r1 = fastq_dir / f"library{r1_suffix}"
                    r2 = fastq_dir / f"library{r2_suffix}"
                    r1.touch()
                    r2.touch()

                    pair = discover_fastq_pair("library", fastq_dir)

                    self.assertEqual(pair.r1, r1.resolve())
                    self.assertEqual(pair.r2, r2.resolve())

    def test_input_id_is_literal_not_a_pattern(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fastq_dir = Path(temporary)
            exact_r1 = fastq_dir / "sample.[1]_R1.fastq.gz"
            exact_r2 = fastq_dir / "sample.[1]_R2.fastq.gz"
            exact_r1.touch()
            exact_r2.touch()
            (fastq_dir / "sample.11_R1.fastq.gz").touch()
            (fastq_dir / "sample.11_R2.fastq.gz").touch()

            pair = discover_fastq_pair("sample.[1]", fastq_dir)

            self.assertEqual(pair.r1, exact_r1.resolve())
            self.assertEqual(pair.r2, exact_r2.resolve())

    def test_does_not_confuse_longer_prefix(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fastq_dir = Path(temporary)
            expected_r1 = fastq_dir / "sample_R1.fastq.gz"
            expected_r2 = fastq_dir / "sample_R2.fastq.gz"
            expected_r1.touch()
            expected_r2.touch()
            (fastq_dir / "sample_extra_R1.fastq.gz").touch()
            (fastq_dir / "sample_extra_R2.fastq.gz").touch()

            pair = discover_fastq_pair("sample", fastq_dir)

            self.assertEqual(pair.r1, expected_r1.resolve())
            self.assertEqual(pair.r2, expected_r2.resolve())

    def test_does_not_search_recursively(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fastq_dir = Path(temporary)
            nested = fastq_dir / "nested"
            nested.mkdir()
            (nested / "sample_R1.fastq.gz").touch()
            (nested / "sample_R2.fastq.gz").touch()

            with self.assertRaisesRegex(FastqNotFoundError, "No paired FASTQ"):
                discover_fastq_pair("sample", fastq_dir)

    def test_path_like_input_id_cannot_escape_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_dir = root / "fastq"
            fastq_dir.mkdir()
            (root / "sample_R1.fastq.gz").touch()
            (root / "sample_R2.fastq.gz").touch()

            with self.assertRaises(FastqNotFoundError):
                discover_fastq_pair("../sample", fastq_dir)

    def test_fails_when_no_pair_exists(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(FastqNotFoundError, "input_id 'missing'"):
                discover_fastq_pair("missing", temporary)

    def test_fails_on_either_paired_end_mismatch(self) -> None:
        cases = (("_R1.fastq.gz", "R2"), ("_R2.fastq.gz", "R1"))
        for present_suffix, missing_mate in cases:
            with self.subTest(missing_mate=missing_mate):
                with tempfile.TemporaryDirectory() as temporary:
                    fastq_dir = Path(temporary)
                    (fastq_dir / f"sample{present_suffix}").touch()

                    with self.assertRaisesRegex(
                        FastqPairMismatchError, f"{missing_mate} is missing"
                    ):
                        discover_fastq_pair("sample", fastq_dir)

    def test_fails_and_lists_all_ambiguous_candidates(self) -> None:
        for ambiguous_mate in ("R1", "R2"):
            with self.subTest(ambiguous_mate=ambiguous_mate):
                with tempfile.TemporaryDirectory() as temporary:
                    fastq_dir = Path(temporary)
                    r1_names = ["sample_R1.fastq.gz"]
                    r2_names = ["sample_R2.fastq.gz"]
                    if ambiguous_mate == "R1":
                        r1_names.append("sample_R1.fq.gz")
                    else:
                        r2_names.append("sample_2.fq.gz")
                    for name in r1_names + r2_names:
                        (fastq_dir / name).touch()

                    with self.assertRaises(AmbiguousFastqError) as raised:
                        discover_fastq_pair("sample", fastq_dir)

                    message = str(raised.exception)
                    for name in r1_names + r2_names:
                        self.assertIn(str((fastq_dir / name).resolve()), message)

    def test_invalid_fastq_directory_is_clear(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            missing_dir = Path(temporary) / "not-there"
            with self.assertRaisesRegex(FastqNotFoundError, "does not exist"):
                discover_fastq_pair("sample", missing_dir)


if __name__ == "__main__":
    unittest.main()
