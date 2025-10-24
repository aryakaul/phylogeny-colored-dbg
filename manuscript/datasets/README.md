# Dataset Acquisition Scripts

This folder contains helper scripts that assemble public genome batches for
pcDBG experiments. The scripts rely on
[`ncbi-genome-download`](https://github.com/kblin/ncbi-genome-download) and only
write to the local filesystem; no data is committed to the repository.

## Workflows

- `download_genomes.py` — downloads high-quality genome cohorts (≥500 assemblies
  per genus by default) into a user-specified output directory and prepares
  matching manifest and tree placeholders.

## Usage

```bash
pip install ncbi-genome-download
python download_genomes.py --output /path/to/datasets
```

The script fetches genomes under `/path/to/datasets/{cohort}/assemblies/`,
selecting the most curated RefSeq assemblies (complete/chromosome level,
reference or representative genomes) until at least 500 sequences are gathered.
It then produces a Snakemake-ready manifest (`{cohort}.txt`) plus stub tree
files that you should replace with real phylogenies.

> Note: `ncbi-genome-download` can retrieve large archives. Use the script's
> `--min-genomes` or `--cohort` flags to scope a run when working on limited
> storage, and remember that metadata resolution requires network access.
