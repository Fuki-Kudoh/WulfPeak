# WulfPeak

WulfPeak uses a compact samplesheet and exact, non-recursive FASTQ discovery.
The v0.2.0 execution foundation supports paired-end (default) and single-end
input validation, full preflight, deterministic dry-run command plans, status
inspection, and output-manifest validation.

This implementation slice intentionally does not execute analysis commands.
Use `wulfpeak run --dry-run`; FastQC through reporting execution is added in
the next implementation PRs.

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

The analysis executables must already be on `PATH`: FastQC, Trim Galore,
Bowtie2, samtools, deepTools `bamCoverage`, MACS3, bedtools, and MultiQC.
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
- `status/pipeline.status` (`null`)
- `status/state.json`

It creates no analysis product and no false `done` status.

## Inspect state and outputs

Status does not require FASTQ, index, or analysis tools:

```console
wulfpeak status --output-dir WulfPeak_out
wulfpeak status --output-dir WulfPeak_out --json
```

`wulfpeak validate-outputs` validates an existing `output_manifest.json`
without rerunning analysis. The complete kind-specific validator set is
finished in the reporting/release PR.

## Development

```console
PYTHONPATH=src PYTHONPYCACHEPREFIX=/tmp/wulfpeak-pycache \
  python -m unittest discover -s tests -v
```
