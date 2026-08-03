"""Command line interface for WulfPeak."""

from __future__ import annotations

import argparse
import sys

from .fastq import FastqDiscoveryError, discover_fastq_pair
from .manifest import ResolvedSample, resolve_sample, write_manifest, write_sample_qc
from .samplesheet import SamplesheetValidationError, load_samplesheet


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="wulfpeak")
    parser.add_argument("--version", action="version", version="%(prog)s 0.1.0")
    subparsers = parser.add_subparsers(dest="command", required=True)

    check = subparsers.add_parser(
        "check", help="validate samples.tsv and resolve every paired FASTQ input"
    )
    check.add_argument(
        "--samplesheet",
        default="samples.tsv",
        help="tab-separated samplesheet (default: samples.tsv)",
    )
    check.add_argument(
        "--fastq-dir",
        default="../fastq",
        help="non-recursive raw FASTQ directory (default: ../fastq)",
    )
    check.add_argument(
        "--output-dir",
        default="WulfPeak_out",
        help="WulfPeak output directory (default: WulfPeak_out)",
    )
    check.set_defaults(func=_run_check)
    return parser


def _run_check(args: argparse.Namespace) -> int:
    samples = load_samplesheet(args.samplesheet)
    resolved_samples: list[ResolvedSample] = []
    discovery_errors: list[str] = []

    for sample in samples:
        try:
            raw_pair = discover_fastq_pair(sample.input_id, args.fastq_dir)
        except FastqDiscoveryError as exc:
            discovery_errors.append(f"sample {sample.sample_id!r}: {exc}")
            continue
        resolved_samples.append(resolve_sample(sample, raw_pair, args.output_dir))

    if discovery_errors:
        raise FastqDiscoveryError(
            "FASTQ discovery failed:\n" + "\n".join(
                f"  - {error}" for error in discovery_errors
            )
        )

    manifest_path = write_manifest(resolved_samples, args.fastq_dir, args.output_dir)
    qc_path = write_sample_qc(resolved_samples, args.output_dir)

    print("WulfPeak check passed (paired-end FASTQ input)")
    for resolved in resolved_samples:
        print(
            f"{resolved.sample.sample_id}:\n"
            f"  R1: {resolved.raw_fastq.r1}\n"
            f"  R2: {resolved.raw_fastq.r2}"
        )
    print(f"Manifest: {manifest_path}")
    print(f"Sample QC: {qc_path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (SamplesheetValidationError, FastqDiscoveryError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
