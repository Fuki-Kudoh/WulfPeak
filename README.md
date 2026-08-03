# WulfPeak

WulfPeak v0.1.0 uses a compact samplesheet and automatically discovers paired
FASTQ files. Single-end reads are not supported in v0.1.0.

## Samplesheet

`samples.tsv` is tab-separated and contains exactly these columns:

| Column | Meaning |
| --- | --- |
| `input_id` | Original sequencing/FASTQ-derived identifier used for raw FASTQ discovery. Must be unique. |
| `sample_id` | Canonical WulfPeak sample name used by all per-sample outputs. Must be unique. |
| `group_id` | Non-empty replicate group used for consensus peaks and pooled BAM/peak generation. |
| `peak_type` | `narrow`, `broad`, or `input`. |
| `control` | Matched input's `sample_id`, or `.`. Input rows must use `.`. |
| `qvalue` | Numeric MACS `callpeak -q` value. Input rows must use `.`. |

See [the example samplesheet](examples/samples.tsv).

The legacy `fastq1` and `fastq2` columns are not accepted. FASTQ paths are
resolved from `input_id` instead.

## Check a run

From the analysis working directory:

```console
wulfpeak check --samplesheet samples.tsv
```

The default raw FASTQ directory is `../fastq`, relative to the current working
directory. Override it when needed:

```console
wulfpeak check --samplesheet samples.tsv --fastq-dir /data/run/fastq
```

`wulfpeak check` validates the entire samplesheet, resolves exactly one R1 and
one R2 for every row, and prints both resolved paths for every sample. On
success it writes:

- `WulfPeak_out/config/manifest.json`
- `WulfPeak_out/qc/sample_qc.tsv`

Both files contain the resolved absolute raw FASTQ paths and canonical trimmed
FASTQ paths. No manifest or sample QC file is written if validation or FASTQ
discovery fails.

## Supported raw FASTQ names

For an `input_id` of `sample`, the accepted pairs are:

```text
sample_R1.fastq.gz       sample_R2.fastq.gz
sample_R1.fq.gz          sample_R2.fq.gz
sample_R1_001.fastq.gz   sample_R2_001.fastq.gz
sample_R1_001.fq.gz      sample_R2_001.fq.gz
sample_1.fastq.gz        sample_2.fastq.gz
sample_1.fq.gz           sample_2.fq.gz
```

Discovery is non-recursive and uses exact literal filenames. If a mate is
missing, no pair exists, or more than one candidate exists for either mate,
the check fails without guessing.

## Trimmed FASTQ contract

Raw names are never used to infer trimmed output names. Each resolved sample in
the manifest declares these canonical downstream paths:

```text
WulfPeak_out/trimmed/{sample_id}_R1_val_1.fq.gz
WulfPeak_out/trimmed/{sample_id}_R2_val_2.fq.gz
```

The trimming step is responsible for moving or linking Trim Galore's produced
files to these paths. Every downstream step must read `trimmed_fastq.r1` and
`trimmed_fastq.r2` from the manifest.

## Development

```console
PYTHONPATH=src python -m unittest discover -s tests -v
```
