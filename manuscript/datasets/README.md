# Dataset Acquisition Scripts

This folder contains helper scripts that assemble public genome batches for
pcDBG experiments. The scripts rely on
[`ncbi-genome-download`](https://github.com/kblin/ncbi-genome-download) and only
write to the local filesystem; no data is committed to the repository.

## Workflows

- `download_genomes.py` — downloads representative collections with varying
  complexity (small, medium, large genome cohorts) into a user-specified
  output directory and prepares matching manifest and tree placeholders.

## Usage

```bash
pip install ncbi-genome-download
python download_genomes.py --output /path/to/datasets
```

The script fetches genomes under `/path/to/datasets/{cohort}/assemblies/` and
produces a Snakemake-ready manifest (`input/{cohort}.txt`) plus stub tree files
that you should replace with real phylogenies.

> Note: `ncbi-genome-download` can retrieve large archives. Adjust the `--limit`
> parameter or edit the script before running on systems with constrained disk
> space.
