from __future__ import annotations

import contextlib
import gzip
import io
import json
import os
import tempfile
import textwrap
import unittest
from pathlib import Path

from wulfpeak.cli import main


INDEX_SUFFIXES = (".1", ".2", ".3", ".4", ".rev.1", ".rev.2")
IMPLEMENTED_TOOLS = ("fastqc", "trim_galore", "bowtie2", "samtools", "bamCoverage")


@contextlib.contextmanager
def environment(**updates: str | None):
    previous = {name: os.environ.get(name) for name in updates}
    for name, value in updates.items():
        if value is None:
            os.environ.pop(name, None)
        else:
            os.environ[name] = value
    try:
        yield
    finally:
        for name, value in previous.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value


FAKE_TOOL = r'''#!/usr/bin/env python3
import gzip
import os
import sys
import zipfile
from pathlib import Path

name = Path(sys.argv[0]).name
args = sys.argv[1:]
if "--version" in args:
    print(f"{name} fake-1.0")
    raise SystemExit(0)

counter = os.environ.get("WULFPEAK_FAKE_COUNTER")
if counter:
    with open(counter, "a", encoding="utf-8") as handle:
        handle.write(name + "\t" + " ".join(args) + "\n")

fail_sample = os.environ.get("WULFPEAK_FAKE_FAIL_SAMPLE")
fail_subcommand = os.environ.get("WULFPEAK_FAKE_FAIL_SUBCOMMAND")
fail = (
    os.environ.get("WULFPEAK_FAKE_FAIL") == name
    and (not fail_sample or fail_sample in " ".join(args))
    and (not fail_subcommand or (args and args[0] == fail_subcommand))
)
invalid = os.environ.get("WULFPEAK_FAKE_INVALID") == name

def option(flag):
    return args[args.index(flag) + 1]

def fastq_name(path, suffix):
    filename = Path(path).name
    for extension in (".fastq.gz", ".fq.gz"):
        if filename.endswith(extension):
            return filename[:-len(extension)] + suffix
    raise RuntimeError(filename)

if name == "fastqc":
    out = Path(option("--outdir"))
    out.mkdir(parents=True, exist_ok=True)
    reads = args[args.index("--outdir") + 2:]
    for index, read in enumerate(reads):
        stem = fastq_name(read, "")
        (out / f"{stem}_fastqc.html").write_text(
            "not html" if invalid else "<!doctype html><html>synthetic</html>\n",
            encoding="utf-8",
        )
        archive = out / f"{stem}_fastqc.zip"
        if invalid:
            archive.write_bytes(b"not-a-zip")
        else:
            with zipfile.ZipFile(archive, "w") as zipped:
                zipped.writestr(f"{stem}_fastqc/fastqc_data.txt", "synthetic\n")
        if fail and index == 0:
            raise SystemExit(7)
    raise SystemExit(0)

if name == "trim_galore":
    out = Path(option("--output_dir"))
    out.mkdir(parents=True, exist_ok=True)
    reads = [value for value in args[args.index("--output_dir") + 2:] if value.endswith(".gz")]
    paired = "--paired" in args
    for index, read in enumerate(reads):
        suffix = ("_val_1.fq.gz" if index == 0 else "_val_2.fq.gz") if paired else "_trimmed.fq.gz"
        destination = out / fastq_name(read, suffix)
        with gzip.open(read, "rb") as source, gzip.open(destination, "wb") as target:
            target.write(source.read())
        (out / (Path(read).name + "_trimming_report.txt")).write_text(
            "synthetic trimming report\n", encoding="utf-8"
        )
    if fail:
        raise SystemExit(8)
    raise SystemExit(0)

if name == "bowtie2":
    if not fail:
        sys.stdout.buffer.write(b"synthetic alignment\n")
    raise SystemExit(9 if fail else 0)

if name == "samtools":
    subcommand = args[0]
    if subcommand == "quickcheck":
        data = Path(args[-1]).read_bytes()
        raise SystemExit(0 if data.startswith(b"BAM") else 1)
    if subcommand == "idxstats":
        bam = Path(args[-1])
        if not bam.is_file() or not bam.read_bytes().startswith(b"BAM"):
            raise SystemExit(1)
        print("chrSynthetic\t1000\t1\t0")
        raise SystemExit(0)
    if fail:
        raise SystemExit(10)
    if subcommand == "view":
        destination = Path(option("-o"))
        sys.stdin.buffer.read()
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(b"BAM\x01synthetic\n")
    elif subcommand == "collate":
        destination = Path(option("-o"))
        sys.stdin.buffer.read()
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(b"BAM\x01synthetic\n")
    elif subcommand == "fixmate":
        data = Path(args[-2]).read_bytes()
        if args[-1] == "-":
            sys.stdout.buffer.write(data)
        else:
            Path(args[-1]).write_bytes(data)
    elif subcommand == "sort":
        data = sys.stdin.buffer.read() if args[-1] == "-" else Path(args[-1]).read_bytes()
        Path(option("-o")).write_bytes(data)
    elif subcommand == "markdup":
        Path(option("-f")).write_text("synthetic markdup metrics\n", encoding="utf-8")
        data = Path(args[-2]).read_bytes()
        if args[-1] == "-":
            sys.stdout.buffer.write(data)
        else:
            Path(args[-1]).write_bytes(data)
    elif subcommand == "index":
        bam = Path(args[-1])
        Path(str(bam) + ".bai").write_bytes(b"BAI\x01synthetic\n")
    elif subcommand in {"flagstat", "stats"}:
        print(f"synthetic {subcommand} output")
    else:
        raise SystemExit(2)
    raise SystemExit(0)

if name == "bamCoverage":
    destination = Path(option("--outFileName"))
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(b"bad!" if invalid else b"\x26\xfc\x8f\x88synthetic\n")
    raise SystemExit(11 if fail else 0)

raise SystemExit(2)
'''


def write_fastq(path: Path, *, records: int = 1) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for index in range(records):
            handle.write(f"@read{index}\nACGT\n+\nIIII\n")


def make_fake_tools(directory: Path) -> None:
    directory.mkdir()
    for name in IMPLEMENTED_TOOLS:
        executable = directory / name
        executable.write_text(textwrap.dedent(FAKE_TOOL), encoding="utf-8")
        executable.chmod(0o755)


def make_index(prefix: Path) -> None:
    for suffix in INDEX_SUFFIXES:
        Path(f"{prefix}{suffix}.bt2").write_bytes(b"synthetic index")


def execution_args(root: Path, *, paired: bool) -> list[str]:
    args = [
        "run",
        "--samplesheet", str(root / "samples.tsv"),
        "--fastq-dir", str(root / "fastq"),
        "--output-dir", str(root / "out"),
        "--assay", "chipseq",
        "--genome-id", "synthetic_reference",
        "--bowtie2-index", str(root / "index" / "reference"),
        "--effective-genome-size", "1000000",
        "--threads", "2",
    ]
    if not paired:
        args.append("--single-end")
    return args


def prepare_run(root: Path, *, paired: bool) -> tuple[list[str], Path]:
    fake_bin = root / "fake_bin"
    make_fake_tools(fake_bin)
    (root / "fastq").mkdir()
    (root / "index").mkdir()
    make_index(root / "index" / "reference")
    (root / "samples.tsv").write_text(
        "input_id\tsample_id\tgroup_id\tpeak_type\tcontrol\tqvalue\n"
        "library_01\tsample_01\tgroup_01\tinput\t.\t.\n",
        encoding="utf-8",
    )
    if paired:
        write_fastq(root / "fastq" / "library_01_R1.fastq.gz")
        write_fastq(root / "fastq" / "library_01_R2.fastq.gz")
    else:
        write_fastq(root / "fastq" / "library_01.fastq.gz")
    return execution_args(root, paired=paired), fake_bin


def prepare_two_sample_run(root: Path) -> tuple[list[str], Path]:
    args, fake_bin = prepare_run(root, paired=False)
    (root / "samples.tsv").write_text(
        "input_id\tsample_id\tgroup_id\tpeak_type\tcontrol\tqvalue\n"
        "library_01\tsample_01\tgroup_01\tinput\t.\t.\n"
        "library_02\tsample_02\tgroup_02\tinput\t.\t.\n",
        encoding="utf-8",
    )
    write_fastq(root / "fastq" / "library_02.fastq.gz")
    return args, fake_bin


def run_with_tools(
    args: list[str], fake_bin: Path, *, counter: Path, fail: str | None = None,
    fail_sample: str | None = None, fail_subcommand: str | None = None,
    invalid: str | None = None,
) -> tuple[int, str]:
    stderr = io.StringIO()
    with environment(
        PATH=f"{fake_bin}{os.pathsep}{os.environ.get('PATH', '')}",
        WULFPEAK_FAKE_COUNTER=str(counter),
        WULFPEAK_FAKE_FAIL=fail,
        WULFPEAK_FAKE_FAIL_SAMPLE=fail_sample,
        WULFPEAK_FAKE_FAIL_SUBCOMMAND=fail_subcommand,
        WULFPEAK_FAKE_INVALID=invalid,
    ), contextlib.redirect_stderr(stderr), contextlib.redirect_stdout(io.StringIO()):
        code = main(args)
    return code, stderr.getvalue()


def tool_counts(counter: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    if not counter.exists():
        return counts
    for line in counter.read_text(encoding="utf-8").splitlines():
        name = line.split("\t", 1)[0]
        counts[name] = counts.get(name, 0) + 1
    return counts


def phase_tool_calls(counter: Path, tool: str, phase: str) -> int:
    if not counter.exists():
        return 0
    return sum(
        1
        for line in counter.read_text(encoding="utf-8").splitlines()
        if line.startswith(f"{tool}\t") and f"/{phase}/" in line
    )


def samtools_analysis_calls(counter: Path) -> int:
    if not counter.exists():
        return 0
    validation_commands = {"quickcheck", "idxstats"}
    return sum(
        1
        for line in counter.read_text(encoding="utf-8").splitlines()
        if line.startswith("samtools\t")
        and line.split("\t", 1)[1].split(" ", 1)[0] not in validation_commands
    )


class ExecutionTests(unittest.TestCase):
    def test_phase_major_fail_fast_resume_continues_with_failed_sample(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args, fake_bin = prepare_two_sample_run(root)
            counter = root / "calls.log"

            code, error = run_with_tools(
                args,
                fake_bin,
                counter=counter,
                fail="trim_galore",
                fail_sample="sample_02",
            )
            self.assertEqual(code, 3)
            self.assertIn("Step trim failed for sample sample_02", error)

            output = root / "out"
            plan = json.loads(
                (output / "config" / "command_plan.json").read_text()
            )
            expected_phases = (
                "fastqc_raw",
                "trim",
                "fastqc_trimmed",
                "align",
                "bam_process",
                "coverage",
            )
            self.assertEqual(
                [
                    (step["phase"], step["scope_id"])
                    for step in plan["steps"]
                ],
                [
                    (phase, sample_id)
                    for phase in expected_phases
                    for sample_id in ("sample_01", "sample_02")
                ],
            )

            state = json.loads((output / "status" / "state.json").read_text())
            self.assertEqual(
                state["steps"]["sample:sample_01:trim"]["status"], "done"
            )
            self.assertEqual(
                state["steps"]["sample:sample_02:trim"]["status"], "failed"
            )
            self.assertFalse(
                any(key.endswith(":fastqc_trimmed") for key in state["steps"])
            )
            status_stdout = io.StringIO()
            with contextlib.redirect_stdout(status_stdout):
                self.assertEqual(
                    main(["status", "--output-dir", str(output)]), 0
                )
            status_text = status_stdout.getvalue()
            self.assertIn("phase: trim", status_text)
            self.assertIn("fastqc_raw:     2/2", status_text)
            self.assertIn("trim:           1/2", status_text)
            self.assertIn("fastqc_trimmed: 0/2", status_text)
            self.assertIn("sample_02  trim", status_text)

            failed_counts = tool_counts(counter)
            self.assertEqual(failed_counts["fastqc"], 2)
            self.assertEqual(failed_counts["trim_galore"], 2)
            failed_raw_fastqc = phase_tool_calls(
                counter, "fastqc", "fastqc_raw"
            )
            code, error = run_with_tools(args, fake_bin, counter=counter)
            self.assertEqual(code, 0, error)
            resumed_counts = tool_counts(counter)
            self.assertEqual(
                phase_tool_calls(counter, "fastqc", "fastqc_raw"),
                failed_raw_fastqc,
            )
            self.assertEqual(
                resumed_counts["trim_galore"], failed_counts["trim_galore"] + 1
            )
            final_state = json.loads(
                (output / "status" / "state.json").read_text()
            )
            self.assertEqual(len(final_state["steps"]), 12)
            self.assertEqual(
                {step["status"] for step in final_state["steps"].values()},
                {"done"},
            )

    def test_bam_scratch_is_cleaned_before_coverage_phase_begins(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args, fake_bin = prepare_two_sample_run(root)
            counter = root / "calls.log"

            code, error = run_with_tools(
                args,
                fake_bin,
                counter=counter,
                fail="samtools",
                fail_sample="sample_02",
                fail_subcommand="fixmate",
            )
            self.assertEqual(code, 3)
            self.assertIn("Step bam_process failed for sample sample_02", error)

            output = root / "out"
            sample1_checkpoint = (
                output
                / "intermediate"
                / "alignment"
                / "sample_01.name-collated.bam"
            )
            sample1_scratch = (
                output
                / "intermediate"
                / ".steps"
                / "samples"
                / "sample_01"
                / "bam_process"
            )
            self.assertFalse(sample1_checkpoint.exists())
            self.assertFalse(sample1_scratch.exists())
            self.assertTrue((output / "bam" / "sample_01.final.bam").is_file())

            sample2_checkpoint = (
                output
                / "intermediate"
                / "alignment"
                / "sample_02.name-collated.bam"
            )
            self.assertTrue(sample2_checkpoint.is_file())
            state = json.loads((output / "status" / "state.json").read_text())
            self.assertEqual(
                state["steps"]["sample:sample_01:bam_process"]["status"],
                "done",
            )
            self.assertEqual(
                state["steps"]["sample:sample_02:bam_process"]["status"],
                "failed",
            )
            self.assertFalse(
                any(key.endswith(":coverage") for key in state["steps"])
            )

    def test_single_end_executes_through_coverage(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args, fake_bin = prepare_run(root, paired=False)
            code, error = run_with_tools(
                [*args, "--stop-after", "coverage"],
                fake_bin,
                counter=root / "calls.log",
            )
            self.assertEqual(code, 0, error)
            output = root / "out"
            metadata = json.loads((output / "metadata" / "run_metadata.json").read_text())
            self.assertFalse(metadata["dry_run"])
            self.assertEqual(metadata["completed_through"], "coverage")
            self.assertEqual((output / "status" / "pipeline.status").read_text().strip(), "done")
            state = json.loads((output / "status" / "state.json").read_text())
            self.assertEqual(len(state["steps"]), 6)
            self.assertEqual({item["status"] for item in state["steps"].values()}, {"done"})
            manifest = json.loads((output / "output_manifest.json").read_text())
            self.assertTrue(manifest["validated"])
            self.assertEqual(manifest["completed_through"], "coverage")
            self.assertEqual({item["phase"] for item in manifest["outputs"]}, {
                "fastqc_raw", "trim", "fastqc_trimmed", "bam_process", "coverage"
            })
            self.assertTrue((output / "bigwig" / "sample_01.normalized.bw").is_file())
            self.assertFalse(
                (
                    output
                    / "intermediate"
                    / "alignment"
                    / "sample_01.name-collated.bam"
                ).exists()
            )
            self.assertFalse(
                (
                    output
                    / "intermediate"
                    / ".steps"
                    / "samples"
                    / "sample_01"
                    / "bam_process"
                ).exists()
            )
            self.assertEqual(
                state["steps"]["sample:sample_01:align"]["output_ids"], []
            )
            selected_path = f"{fake_bin}{os.pathsep}{os.environ.get('PATH', '')}"
            with environment(PATH=selected_path), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(
                    main(["validate-outputs", "--output-dir", str(output)]), 0
                )

    def test_paired_end_executes_and_reruns_nonempty_fastqc_directories(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args, fake_bin = prepare_run(root, paired=True)
            counter = root / "calls.log"
            code, error = run_with_tools(args, fake_bin, counter=counter)
            self.assertEqual(code, 0, error)
            output = root / "out"
            raw_dir = output / "qc" / "fastqc" / "raw" / "sample_01"
            trimmed_dir = output / "qc" / "fastqc" / "trimmed" / "sample_01"
            for directory in (raw_dir, trimmed_dir):
                (directory / "sentinel.txt").write_text("keep\n", encoding="utf-8")
            code, error = run_with_tools(
                [*args, "--force-from", "fastqc_raw"], fake_bin, counter=counter
            )
            self.assertEqual(code, 0, error)
            self.assertEqual((raw_dir / "sentinel.txt").read_text(), "keep\n")
            self.assertEqual((trimmed_dir / "sentinel.txt").read_text(), "keep\n")
            manifest = json.loads((output / "output_manifest.json").read_text())
            fastqc_outputs = [item for item in manifest["outputs"] if item["phase"].startswith("fastqc")]
            self.assertEqual(len(fastqc_outputs), 8)
            self.assertTrue(all(item["validated"] for item in fastqc_outputs))

    def test_invalid_fastqc_set_promotes_none_and_marks_failure(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args, fake_bin = prepare_run(root, paired=False)
            canonical = root / "out" / "qc" / "fastqc" / "raw" / "sample_01"
            canonical.mkdir(parents=True)
            html = canonical / "library_01_fastqc.html"
            archive = canonical / "library_01_fastqc.zip"
            html.write_text("old html\n", encoding="utf-8")
            archive.write_bytes(b"old zip")
            code, error = run_with_tools(
                args, fake_bin, counter=root / "calls.log", invalid="fastqc"
            )
            self.assertEqual(code, 3)
            self.assertIn("Temporary artifact validation failed", error)
            self.assertEqual(html.read_text(), "old html\n")
            self.assertEqual(archive.read_bytes(), b"old zip")
            output = root / "out"
            self.assertEqual((output / "status" / "pipeline.status").read_text().strip(), "failed")
            state = json.loads((output / "status" / "state.json").read_text())
            self.assertEqual(state["steps"]["sample:sample_01:fastqc_raw"]["status"], "failed")

    def test_command_failure_never_promotes_partial_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args, fake_bin = prepare_run(root, paired=False)
            canonical = root / "out" / "qc" / "fastqc" / "raw" / "sample_01"
            canonical.mkdir(parents=True)
            html = canonical / "library_01_fastqc.html"
            html.write_text("previous\n", encoding="utf-8")
            code, _error = run_with_tools(
                args, fake_bin, counter=root / "calls.log", fail="fastqc"
            )
            self.assertEqual(code, 3)
            self.assertEqual(html.read_text(), "previous\n")
            state = json.loads(
                (root / "out" / "status" / "state.json").read_text()
            )
            self.assertEqual(
                state["steps"]["sample:sample_01:fastqc_raw"]["exit_code"], 7
            )

    def test_resume_skips_all_commands_and_input_change_rebuilds(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args, fake_bin = prepare_run(root, paired=False)
            counter = root / "calls.log"
            self.assertEqual(run_with_tools(args, fake_bin, counter=counter)[0], 0)
            first = tool_counts(counter)
            first_samtools_analysis = samtools_analysis_calls(counter)
            self.assertEqual(run_with_tools(args, fake_bin, counter=counter)[0], 0)
            resumed = tool_counts(counter)
            for tool in ("fastqc", "trim_galore", "bowtie2", "bamCoverage"):
                self.assertEqual(resumed[tool], first[tool])
            self.assertEqual(
                samtools_analysis_calls(counter), first_samtools_analysis
            )
            write_fastq(root / "fastq" / "library_01.fastq.gz", records=2)
            self.assertEqual(run_with_tools(args, fake_bin, counter=counter)[0], 0)
            rebuilt = tool_counts(counter)
            self.assertGreater(rebuilt["fastqc"], resumed["fastqc"])
            self.assertGreater(rebuilt["bamCoverage"], resumed["bamCoverage"])

    def test_force_from_and_no_resume_select_expected_steps(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args, fake_bin = prepare_run(root, paired=False)
            counter = root / "calls.log"
            self.assertEqual(run_with_tools(args, fake_bin, counter=counter)[0], 0)
            first = tool_counts(counter)
            self.assertEqual(
                run_with_tools([*args, "--force-from", "align"], fake_bin, counter=counter)[0],
                0,
            )
            forced = tool_counts(counter)
            self.assertEqual(forced["fastqc"], first["fastqc"])
            self.assertEqual(forced["trim_galore"], first["trim_galore"])
            self.assertGreater(forced["bowtie2"], first["bowtie2"])
            self.assertGreater(forced["bamCoverage"], first["bamCoverage"])
            self.assertEqual(
                run_with_tools(
                    [*args, "--force-from", "bam_process"],
                    fake_bin,
                    counter=counter,
                )[0],
                0,
            )
            forced_bam = tool_counts(counter)
            self.assertEqual(forced_bam["fastqc"], forced["fastqc"])
            self.assertEqual(forced_bam["trim_galore"], forced["trim_galore"])
            self.assertGreater(forced_bam["bowtie2"], forced["bowtie2"])
            self.assertGreater(forced_bam["bamCoverage"], forced["bamCoverage"])
            self.assertEqual(
                run_with_tools([*args, "--no-resume"], fake_bin, counter=counter)[0], 0
            )
            no_resume = tool_counts(counter)
            self.assertGreater(no_resume["fastqc"], forced_bam["fastqc"])
            self.assertGreater(no_resume["trim_galore"], forced_bam["trim_galore"])

    def test_missing_final_bam_rehydrates_retired_alignment(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args, fake_bin = prepare_run(root, paired=False)
            counter = root / "calls.log"
            self.assertEqual(run_with_tools(args, fake_bin, counter=counter)[0], 0)
            initial = tool_counts(counter)
            final_bam = root / "out" / "bam" / "sample_01.final.bam"
            final_bam.unlink()

            code, error = run_with_tools(args, fake_bin, counter=counter)
            self.assertEqual(code, 0, error)
            rebuilt = tool_counts(counter)
            self.assertGreater(rebuilt["bowtie2"], initial["bowtie2"])
            self.assertGreater(rebuilt["bamCoverage"], initial["bamCoverage"])
            self.assertTrue(final_bam.is_file())
            self.assertFalse(
                (
                    root
                    / "out"
                    / "intermediate"
                    / "alignment"
                    / "sample_01.name-collated.bam"
                ).exists()
            )

    def test_invalid_forced_coverage_preserves_previous_bigwig(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args, fake_bin = prepare_run(root, paired=False)
            counter = root / "calls.log"
            self.assertEqual(run_with_tools(args, fake_bin, counter=counter)[0], 0)
            initial = tool_counts(counter)
            bigwig = root / "out" / "bigwig" / "sample_01.normalized.bw"
            previous = bigwig.read_bytes()
            code, _error = run_with_tools(
                [*args, "--force-from", "coverage"],
                fake_bin,
                counter=counter,
                invalid="bamCoverage",
            )
            self.assertEqual(code, 3)
            self.assertEqual(bigwig.read_bytes(), previous)
            self.assertEqual(tool_counts(counter)["bowtie2"], initial["bowtie2"])

    def test_bam_process_cleans_scratch_before_coverage_failure(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args, fake_bin = prepare_run(root, paired=False)
            counter = root / "calls.log"
            code, _error = run_with_tools(
                args,
                fake_bin,
                counter=counter,
                invalid="bamCoverage",
            )
            self.assertEqual(code, 3)
            output = root / "out"
            checkpoint = (
                output
                / "intermediate"
                / "alignment"
                / "sample_01.name-collated.bam"
            )
            scratch = (
                output
                / "intermediate"
                / ".steps"
                / "samples"
                / "sample_01"
                / "bam_process"
            )
            self.assertFalse(checkpoint.exists())
            self.assertFalse(scratch.exists())
            self.assertTrue((output / "bam" / "sample_01.final.bam").is_file())
            before_resume = tool_counts(counter)
            before_samtools = samtools_analysis_calls(counter)

            self.assertEqual(run_with_tools(args, fake_bin, counter=counter)[0], 0)
            after_resume = tool_counts(counter)
            self.assertEqual(after_resume["bowtie2"], before_resume["bowtie2"])
            self.assertEqual(samtools_analysis_calls(counter), before_samtools)

    def test_keep_intermediates_retains_collated_checkpoint_and_sort_scratch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args, fake_bin = prepare_run(root, paired=False)
            counter = root / "calls.log"
            self.assertEqual(
                run_with_tools(
                    [*args, "--keep-intermediates"], fake_bin, counter=counter
                )[0],
                0,
            )
            output = root / "out"
            checkpoint = (
                output
                / "intermediate"
                / "alignment"
                / "sample_01.name-collated.bam"
            )
            self.assertTrue(checkpoint.is_file())
            scratch = (
                output
                / "intermediate"
                / ".steps"
                / "samples"
                / "sample_01"
                / "bam_process"
            )
            self.assertEqual(
                {
                    path.name
                    for path in scratch.iterdir()
                    if path.suffix == ".bam"
                },
                {"coordinate.bam"},
            )
            manifest = json.loads((output / "output_manifest.json").read_text())
            align_outputs = [
                entry for entry in manifest["outputs"] if entry["phase"] == "align"
            ]
            self.assertEqual(len(align_outputs), 1)
            self.assertEqual(
                Path(align_outputs[0]["path"]).resolve(), checkpoint.resolve()
            )

            analysis_calls = {
                tool: tool_counts(counter)[tool]
                for tool in ("fastqc", "trim_galore", "bowtie2", "bamCoverage")
            }
            samtools_calls = samtools_analysis_calls(counter)
            self.assertEqual(
                run_with_tools(
                    [*args, "--keep-intermediates"], fake_bin, counter=counter
                )[0],
                0,
            )
            self.assertEqual(
                {
                    tool: tool_counts(counter)[tool]
                    for tool in ("fastqc", "trim_galore", "bowtie2", "bamCoverage")
                },
                analysis_calls,
            )
            self.assertEqual(samtools_analysis_calls(counter), samtools_calls)

            self.assertEqual(run_with_tools(args, fake_bin, counter=counter)[0], 0)
            self.assertEqual(
                {
                    tool: tool_counts(counter)[tool]
                    for tool in ("fastqc", "trim_galore", "bowtie2", "bamCoverage")
                },
                analysis_calls,
            )
            self.assertFalse(checkpoint.exists())
            self.assertFalse(scratch.exists())
            manifest = json.loads((output / "output_manifest.json").read_text())
            self.assertNotIn(
                "align", {entry["phase"] for entry in manifest["outputs"]}
            )


if __name__ == "__main__":
    unittest.main()
