# Compression Metrics Analysis

This analysis summarizes the lossy color compression that occurs when moving
from the original Cuttlefish-colored compacted DBG to the phylogeny-colored
DBG produced by pcDBG. It consumes the metrics emitted by the workflow under
`intermediate/metrics/` and generates manuscript-ready tables summarizing the
number of unique color sets before and after parsimony.

## Prerequisites

1. Run the pcDBG workflow (e.g., `make all`) with the desired batches and
   parsimony configuration.
2. Ensure the compression metrics TSVs exist:
   ```
   intermediate/metrics/{batch}_compression_{mode_label}_k{k}.tsv
   ```
   Each file contains tab-separated `metric`, `value`, and `description` fields.

## Running the analysis

Use `analyze_metrics.py` to aggregate one or more metric files:

```bash
python analyze_metrics.py \
  --metrics intermediate/metrics/mtbc_compression_fitch_k31.tsv \
  --metrics intermediate/metrics/mtbc_compression_sankoff-auto_k31.tsv \
  --output results/compression_summary.tsv
```

The script accepts multiple `--metrics` arguments and writes a tabulated summary
of the unique color counts, absolute differences, and ratios for each file. If
`--output` is omitted the table prints to stdout.

## Outputs

- Summary table listing the batch, parsimony label, k-mer length, and the four
  compression metrics extracted from each TSV.
- Optional plots or manuscript tables can be generated downstream using the
  summary as input (not included here to avoid large binary artifacts).

## Reproducibility notes

- The analysis script performs pure aggregation and does not require access to
  raw assemblies or phylogenies; reproducing the TSVs requires rerunning the
  Snakemake workflow with the same configuration.
- If additional diagnostics (e.g., sample-level color distributions) are
  needed, extend this directory with separate scripts and document any
  additional required workflow outputs.
