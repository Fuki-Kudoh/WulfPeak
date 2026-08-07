# WulfPeak

WulfPeak uses a compact samplesheet and exact, non-recursive FASTQ discovery.
The v0.2.1 implementation supports paired-end (default) and single-end input
validation, deterministic command plans, status inspection, conservative
resume, real per-sample execution through normalized BigWig coverage, and a
pipeline-scoped MultiQC report after all samples complete coverage.

## Samplesheet

`samples.tsv` is tab-separated and contains exactly these columns:

| Column | Meaning |
| --- | --- |
| `input_id` | Raw FASTQ discovery identifier. Unique within the run. |
| `sample_id` | Path-safe canonical output name. Unique within the run. |
| `group_id` | Path-safe replicate group identifier. |
| `peak_type` | `narrow`, `broad`, or `input`. |
| `control` | Matched input `sample_id`, or `.`. Input rows use `.`. |
| `qvalue` | Finite `0 < qvalue <= 1` for treatment; `.` for input. |

See [the paired-end example](examples/samples.tsv) and
[single-end example](examples/samples-se.tsv). The samplesheet schema is the
same for both layouts.

Treatment replicates in one `group_id` must have the same peak type,
numerically equivalent q-value, and control. The legacy `fastq1` and `fastq2`
columns are not accepted.

## Check a run

Paired-end is the default:

```console
wulfpeak check --samplesheet samples.tsv --fastq-dir ../fastq
```

Select single-end with either alias:

```console
wulfpeak check --samplesheet samples.tsv --single-end
wulfpeak check --samplesheet samples.tsv --SE
```

On success, `check` writes `WulfPeak_out/config/manifest.json` and
`WulfPeak_out/qc/sample_qc.tsv`. It writes neither until every row passes.
Single-end manifests use JSON `null` for R2 and the TSV R2 columns are empty.

## Supported raw FASTQ names

For an `input_id` of `sample`, accepted paired-end pairs are:

```text
sample_R1.fastq.gz       sample_R2.fastq.gz
sample_R1.fq.gz          sample_R2.fq.gz
sample_R1_001.fastq.gz   sample_R2_001.fastq.gz
sample_R1_001.fq.gz      sample_R2_001.fq.gz
sample_1.fastq.gz        sample_2.fastq.gz
sample_1.fq.gz           sample_2.fq.gz
```

Accepted single-end names are:

```text
sample.fastq.gz
sample.fq.gz
sample_R1.fastq.gz
sample_R1.fq.gz
sample_R1_001.fastq.gz
sample_R1_001.fq.gz
sample_1.fastq.gz
sample_1.fq.gz
```

Discovery is exact and non-recursive. Missing, orphaned, or multiple candidates
fail without guessing. R2 is an error in single-end mode, never silently
ignored.

## Canonical trimmed FASTQ paths

Paired-end:

```text
WulfPeak_out/trimmed/{sample_id}_R1_val_1.fq.gz
WulfPeak_out/trimmed/{sample_id}_R2_val_2.fq.gz
```

Single-end:

```text
WulfPeak_out/trimmed/{sample_id}_trimmed.fq.gz
```

Downstream steps consume these manifest paths rather than reconstructing names.

## Validate and plan a run

The executables required for the default real execution through MultiQC must
already be on `PATH`: FastQC, Trim Galore, Bowtie2, samtools, deepTools
`bamCoverage`, and MultiQC. A run with `--stop-after coverage` does not require
MultiQC. A complete dry-run plan additionally requires MACS3 and bedtools.
WulfPeak does not install tools or load environment modules.

```console
wulfpeak run \
  --samplesheet samples.tsv \
  --fastq-dir ../fastq \
  --output-dir WulfPeak_out \
  --assay chipseq \
  --genome-id synthetic_reference \
  --bowtie2-index /path/to/index/prefix \
  --effective-genome-size 1000000 \
  --threads 8 \
  --dry-run
```

Dry-run validates the full samplesheet, FASTQ files, numeric options, complete
Bowtie2 index family, optional blacklist, output directory, and required tools.
It writes:

- `config/manifest.json`
- `config/command_plan.json`
- `metadata/run_metadata.json`
- `status/pipeline.status` (`null` for a new output directory)
- `status/state.json` (initialized only for a new output directory)

It creates no analysis product and no false `done` status. Re-running dry-run
against an existing output directory preserves all completed and failed state.

Each planned step has ordered `actions` and explicit `artifacts`. Commands that
produce data on stdout declare an atomic `stdout_path`. Every artifact records
a distinct step-temporary path, validator, canonical path, and promotion
contract. Files use `validate_then_atomic_file_replace`. The complete MultiQC
directory uses a validated, same-parent staged directory replacement, avoiding
cross-filesystem rename assumptions and preventing partial reports from being
promoted. The validated temporary directory is removed after successful
promotion. FastQC and MACS3 outputs are declared as individual files, so a
forced rerun never replaces a non-empty parent directory. The consensus plan
normalizes each replicate to merged BED3, records strict-majority support,
merges qualifying segments, and declares both the canonical BED3 and support
TSV outputs.

## Execute through MultiQC

Omitting `--dry-run` executes all applicable samples in phase-major order, with
a barrier between phases:

```text
fastqc_raw:      all samples
trim:            all samples
fastqc_trimmed:  all samples
align:           all samples
bam_process:     all samples
coverage:        all samples
multiqc:         pipeline
```

A later phase never starts until the current phase has completed for every
applicable sample. Execution remains fail-fast: one sample failure stops the
run, and resume skips reusable completed samples in that phase before
continuing from failed or incomplete work.

`--threads` is an approximate per-sample CPU budget. Concurrent pipelines use
a centralized allocation policy instead of assigning that full value to every
process: alignment divides the budget between Bowtie2 and the samtools collate
main/worker threads, while BAM processing reserves one main thread each for
fixmate and sort and assigns the remainder to sort workers. A budget of one
uses no additional samtools workers, although a streaming pipeline still
requires its two main processes. The resolved allocation is recorded in
`command_plan.json` as `thread_allocation`.

`--stop-after` defaults to `multiqc` for a real run. Use `--stop-after coverage`
to stop before consolidated reporting, or `--stop-after multiqc` explicitly for
the default v0.2.1 boundary. A later stop such as `peak` fails before preflight.
MACS3 peak calling, pooled BAMs, pooled and consensus peaks, and the final
WulfPeak `report` phase remain plan-only and are never executed by this release.
The final `report` phase is distinct from MultiQC and remains reserved for a
future WulfPeak run summary.

Each step writes into its owned temporary directory. All temporary artifacts
must pass their declared validators before any artifact in the step is
atomically promoted. FastQC HTML and ZIP files are promoted individually, so
reruns retain their existing non-empty parent directories.

Resume is enabled by default. A step is reused only when its prior state is
`done`, its signature still matches the plan, tool provenance, options, and
direct input fingerprints, every retained canonical output validates, and its
upstream chain was also reusable. Use `--force-from PHASE` to rebuild that
phase and everything after it, or `--no-resume` to rebuild all implemented
phases. MultiQC is signed from every sample's coverage chain: a failed MultiQC
step resumes without rebuilding BAMs or BigWigs, while rebuilt upstream
coverage invalidates MultiQC normally.

Alignment streams directly from Bowtie2 into a name-collated BAM, which is the
resumable alignment checkpoint. BAM processing streams `fixmate` into `sort`
and `markdup` into the final filtering step, avoiding unsorted, fixmate, and
marked BAM materializations. By default, WulfPeak removes the collated
checkpoint and coordinate-sort scratch as soon as that sample completes
`bam_process`; coverage failures therefore do not retain cohort-wide BAM
scratch. Resume reuses the downstream validated final BAM or regenerates the
alignment checkpoint when BAM processing must rerun. Pass `--keep-intermediates`
to retain the collated checkpoint and coordinate-sort scratch; the checkpoint
is then included in `output_manifest.json`.

Normalized BigWig generation remains the responsibility of `coverage`; v0.2.1
does not change that processing. It adds consolidated QC reporting from the
FastQC, Trim Galore, and samtools outputs already present under the canonical
`qc/` directory. MultiQC scans only this canonical QC tree, excluding step
temporary files and reports from previous runs. Original Trim Galore text and
v2 JSON reports are preserved under `qc/trimming/<sample>/` so MultiQC never
depends on files under `intermediate/.steps/`. A successful default run
preserves the standard MultiQC outputs:

```text
multiqc/
├── multiqc_report.html
└── multiqc_data/
```

Other standard MultiQC companion files are preserved as well. A successful
default real run records `dry_run: false` and `completed_through: multiqc` in
run metadata and `output_manifest.json`, sets all six per-sample step states to
`done`, and records the single pipeline-scoped MultiQC step as `done`.

## Inspect state and outputs

Status does not require FASTQ, index, or analysis tools:

```console
wulfpeak status --output-dir WulfPeak_out
wulfpeak status --output-dir WulfPeak_out --json
```

The human-readable status includes the current phase, completed/total counts
for each planned phase, and details for running or failed work. Per-sample
phases count against the number of samples; pipeline-scoped MultiQC is shown as
`0/1` or `1/1` and its detail is identified as `pipeline`, not as a sample.

`wulfpeak validate-outputs` validates an existing `output_manifest.json`
without rerunning analysis. Validation checks gzip FASTQ structure, FastQC
HTML and ZIP integrity, final BAMs with `samtools quickcheck`, BAM indexes,
non-empty QC text, Trim Galore source reports, BigWig magic bytes, and the
MultiQC report/data directory. MultiQC validation also requires a non-empty
`multiqc_data/multiqc_sources.txt`, confirming that at least one source was
parsed. Because BAM validation is explicit,
`samtools` must be available when a manifest contains BAM artifacts.

## Development

```console
PYTHONPATH=src PYTHONPYCACHEPREFIX=/tmp/wulfpeak-pycache \
  python -m unittest discover -s tests -v
```
