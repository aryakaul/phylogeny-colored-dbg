# Manuscript Analyses

This directory tracks reproducible analysis recipes for the pcDBG manuscript.
Each subfolder bundles short instructions with scripts that consume workflow
outputs; rerun Snakemake with the needed configuration before using them.

Available analyses:

- [`compression_metrics/`](compression_metrics/) — contrasts color compression
  between the original Cuttlefish compacted DBG and the phylogeny-colored DBG.

General guidelines:

1. Run `make all` (or the specific Snakemake targets) so required outputs
   exist before launching an analysis.
2. Keep large results out of the repository or add them to `.gitignore`.
3. Make scripts parameterized so readers can reuse them on their own data.
