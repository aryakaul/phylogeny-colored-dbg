# Dataset Acquisition Scripts

This folder contains helper scripts that assemble public genome batches for
pcDBG experiments. The scripts rely on
[`ncbi-genome-download`](https://github.com/kblin/ncbi-genome-download) and only
write to the local filesystem; no data is committed to the repository.

## Workflows

- `get_sets.py` — samples three matched cohorts (SetL/SetM/SetH) at increasing
  diversity levels for compression metrics, caching assembly metadata, FASTA
  downloads, and Mash sketches.

## Usage

```bash
python get_sets.py --n 50 --threads 8
```

`get_sets.py` builds three genome sets of increasing diversity (SetL/SetM/SetH),
caching NCBI assembly metadata, FASTA downloads and Mash sketches between runs.
See `python get_sets.py --help` for the full set of options, including plasmid
handling, Mash sketch parameters, and the distance thresholds used to define
SetH.

> Note: `ncbi-genome-download` can retrieve large archives. Use the script's
> `--min-genomes` or `--cohort` flags to scope a run when working on limited
> storage, and remember that metadata resolution requires network access.

## Sampling Low/Medium/High Diversity Sets

`get_sets.py` constructs three size-matched cohorts (default `N=500`) under
`work/`, `sets/`, and `figs/`:

- **Set L (low diversity)** — chooses the RefSeq species with the deepest pool
  of latest complete/chromosome assemblies, then auto-tunes a Mash distance
  threshold (starting at 0.0015) to extract `N` genomes from the densest
  within-species component. Length decile quotas are recorded for downstream
  matching.
- **Set M (medium diversity)** — runs staged farthest-first selection inside the
  same species with a Mash floor of 0.005 and cap of 0.06. If quotas prevent a
  full set, the floor relaxes toward 0.002 before a final capped fill ensures
  `N` genomes while keeping length deciles aligned with Set L.
- **Set H (high diversity)** — curates one representative/reference assembly per
  genus (optionally restricted to a phylum via `--h-phylum` and `--taxdump`),
  then applies farthest-first with a 0.10 floor and 0.95 cap, padding with
  moderate-distance genomes if needed.

All stages cache `ncbi-genome-download` outputs under `work/candidates_*` and
reuse Mash sketches unless `--force-*` flags are supplied. Outputs include
symlinked FASTAs per set, `sets/summary.tsv`, taxonomy tables
(`sets/df_tax.tsv/.parquet` when taxdump data is provided), and seaborn plots in
`figs/`.

### Example

```bash
pip install ncbi-genome-download mash seaborn pandas numpy
python get_sets.py --n 500 --threads 8
```

Use `--length-match/--no-length-match` to toggle length quotas, `--plasmids` to
drop plasmid contigs, and `--force-metadata/--force-download/--force-sketch` to
refresh caches when new inputs are desired. All downloads remain outside the
repository; keep an eye on disk usage before scaling `N`.
