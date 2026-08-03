"""Command line interface for WulfPeak."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from . import __version__
from .command import CommandExecutionError
from .config import PHASES, ConfigurationError, RunConfig
from .fastq import FastqDiscoveryError, discover_fastq_input
from .manifest import ResolvedSample, resolve_sample, write_manifest, write_sample_qc
from .models import Assay, ReadLayout
from .runner import prepare_dry_run
from .samplesheet import SamplesheetValidationError, load_samplesheet
from .status import LockConflictError, RunLock, read_status
from .validators import OutputValidationError, validate_output_manifest


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def _nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be zero or greater")
    return parsed


def _cutoff(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed) or not 0 < parsed <= 1:
        raise argparse.ArgumentTypeError("must be finite and satisfy 0 < value <= 1")
    return parsed


def _add_layout(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--single-end",
        "--SE",
        dest="single_end",
        action="store_true",
        help="use single-end input (paired-end is the default)",
    )


def _add_common_inputs(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--samplesheet", default="samples.tsv")
    parser.add_argument("--fastq-dir", default="../fastq")
    parser.add_argument("--output-dir", default="WulfPeak_out")
    _add_layout(parser)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="wulfpeak")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    check = subparsers.add_parser(
        "check", help="validate samplesheet and resolve FASTQ inputs"
    )
    _add_common_inputs(check)
    check.set_defaults(func=_run_check)

    run = subparsers.add_parser(
        "run", help="preflight and run (execution foundation: dry-run)"
    )
    _add_common_inputs(run)
    run.add_argument("--assay", required=True, choices=[item.value for item in Assay])
    run.add_argument("--genome-id", required=True)
    run.add_argument("--bowtie2-index", required=True)
    run.add_argument("--effective-genome-size", required=True, type=_positive_int)
    run.add_argument("--threads", required=True, type=_positive_int)
    run.add_argument("--min-mapq", type=_nonnegative_int, default=30)
    run.add_argument("--duplicate-policy", choices=("remove", "keep"), default="remove")
    run.add_argument("--allow-dovetail", action="store_true")
    run.add_argument("--normalization", choices=("CPM", "RPGC", "None"), default="CPM")
    run.add_argument("--bin-size", type=_positive_int, default=10)
    run.add_argument("--blacklist")
    run.add_argument("--broad-cutoff", type=_cutoff, default=0.1)
    run.add_argument("--keep-intermediates", action="store_true")
    run.add_argument("--dry-run", action="store_true")
    run.add_argument("--no-resume", action="store_true")
    run.add_argument("--force-from", choices=PHASES)
    run.set_defaults(func=_run_pipeline)

    status = subparsers.add_parser(
        "status", help="read pipeline status without analysis tools"
    )
    status.add_argument("--output-dir", default="WulfPeak_out")
    status.add_argument("--json", action="store_true")
    status.set_defaults(func=_run_status)

    validate = subparsers.add_parser(
        "validate-outputs", help="validate declared outputs only"
    )
    validate.add_argument("--output-dir", default="WulfPeak_out")
    validate.add_argument("--json", action="store_true")
    validate.set_defaults(func=_run_validate_outputs)
    return parser


def _layout(args: argparse.Namespace) -> ReadLayout:
    return ReadLayout.SINGLE_END if args.single_end else ReadLayout.PAIRED_END


def _run_check(args: argparse.Namespace) -> int:
    samples = load_samplesheet(args.samplesheet)
    layout = _layout(args)
    resolved_samples: list[ResolvedSample] = []
    discovery_errors: list[str] = []
    for sample in samples:
        try:
            raw = discover_fastq_input(sample.input_id, args.fastq_dir, layout)
        except FastqDiscoveryError as exc:
            discovery_errors.append(f"sample {sample.sample_id!r}: {exc}")
            continue
        resolved_samples.append(resolve_sample(sample, raw, args.output_dir, layout))
    if discovery_errors:
        raise FastqDiscoveryError(
            "FASTQ discovery failed:\n" + "\n".join(
                f"  - {error}" for error in discovery_errors
            )
        )
    manifest_path = write_manifest(
        resolved_samples,
        args.fastq_dir,
        args.output_dir,
        layout,
        args.samplesheet,
    )
    qc_path = write_sample_qc(resolved_samples, args.output_dir)
    print(f"WulfPeak check passed ({layout.value})")
    for resolved in resolved_samples:
        print(f"{resolved.sample.sample_id}:\n  R1: {resolved.raw_fastq.r1}")
        if resolved.raw_fastq.r2 is not None:
            print(f"  R2: {resolved.raw_fastq.r2}")
    print(f"Manifest: {manifest_path}\nSample QC: {qc_path}")
    return 0


def _run_pipeline(args: argparse.Namespace) -> int:
    config = RunConfig(
        samplesheet=Path(args.samplesheet).expanduser().resolve(),
        fastq_dir=Path(args.fastq_dir).expanduser().resolve(),
        output_dir=Path(args.output_dir).expanduser().resolve(),
        assay=Assay(args.assay),
        genome_id=args.genome_id,
        bowtie2_index=Path(args.bowtie2_index).expanduser().resolve(),
        effective_genome_size=args.effective_genome_size,
        threads=args.threads,
        read_layout=_layout(args),
        min_mapq=args.min_mapq,
        duplicate_policy=args.duplicate_policy,
        allow_dovetail=args.allow_dovetail,
        normalization=args.normalization,
        bin_size=args.bin_size,
        blacklist=(
            Path(args.blacklist).expanduser().resolve() if args.blacklist else None
        ),
        broad_cutoff=args.broad_cutoff,
        keep_intermediates=args.keep_intermediates,
        dry_run=args.dry_run,
        resume=not args.no_resume,
        force_from=args.force_from,
    )
    manifest, plan, metadata = prepare_dry_run(config, args._argv)
    print("WulfPeak dry-run preflight passed; no analysis commands were executed")
    print(f"Manifest: {manifest}\nCommand plan: {plan}\nRun metadata: {metadata}")
    return 0


def _run_status(args: argparse.Namespace) -> int:
    payload = read_status(args.output_dir)
    if args.json:
        print(json.dumps(payload, sort_keys=True))
    else:
        print(f"pipeline: {payload['pipeline_display']}")
        for key, state in sorted(payload.get("steps", {}).items()):
            print(f"{key}: {state.get('display_status', state.get('status', 'null'))}")
    return 0


def _run_validate_outputs(args: argparse.Namespace) -> int:
    with RunLock(args.output_dir, ["wulfpeak", "validate-outputs"]):
        path, valid, payload = validate_output_manifest(args.output_dir)
    if args.json:
        print(json.dumps(payload, sort_keys=True))
    else:
        print(f"Output manifest: {path}\nvalid: {str(valid).lower()}")
    return 0 if valid else 4


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    args._argv = ["wulfpeak", *(argv if argv is not None else sys.argv[1:])]
    try:
        return args.func(args)
    except LockConflictError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 5
    except CommandExecutionError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 3
    except OutputValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 4
    except (
        SamplesheetValidationError,
        FastqDiscoveryError,
        ConfigurationError,
        ValueError,
        OSError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
