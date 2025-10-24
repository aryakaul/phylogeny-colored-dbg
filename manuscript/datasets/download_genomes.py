#!/usr/bin/env python3
"""Download representative genome cohorts for pcDBG experiments.

The script wraps `ncbi-genome-download` to pull a small number of genomes for
species with varying genomic complexity. It writes manifests that can be fed
directly into the pcDBG workflow and creates placeholder tree files for manual
curation.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from shutil import which
from typing import Iterable, List, Optional


@dataclass(frozen=True)
class Cohort:
    identifier: str
    group: str
    species: str
    default_limit: int


COHORTS: List[Cohort] = [
    Cohort(
        identifier="listeria_low",
        group="bacteria",
        species="Listeria monocytogenes",
        default_limit=10,
    ),
    Cohort(
        identifier="mtbc_medium",
        group="bacteria",
        species="Mycobacterium tuberculosis",
        default_limit=20,
    ),
    Cohort(
        identifier="ecoli_high",
        group="bacteria",
        species="Escherichia coli",
        default_limit=40,
    ),
]


def parse_args(argv: Optional[Iterable[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download genome cohorts using ncbi-genome-download."
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Base directory where cohorts, manifests, and trees will be written.",
    )
    parser.add_argument(
        "--formats",
        default="fasta",
        help="Comma-separated list of formats (passed to ncbi-genome-download).",
    )
    parser.add_argument(
        "--section",
        default="refseq",
        help="NCBI section to download from (refseq or genbank).",
    )
    parser.add_argument(
        "--limit",
        type=int,
        help="Override the number of genomes to fetch for every cohort.",
    )
    parser.add_argument(
        "--cohort",
        action="append",
        choices=[c.identifier for c in COHORTS],
        help="Download only the specified cohort(s). Can be used multiple times.",
    )
    return parser.parse_args(argv)


def ensure_tool_present() -> None:
    if which("ncbi-genome-download") is None:
        raise RuntimeError(
            "ncbi-genome-download is not available on PATH. "
            "Install it via `pip install ncbi-genome-download`."
        )


def run_download(
    cohort: Cohort,
    out_dir: Path,
    formats: str,
    section: str,
    limit_override: Optional[int],
) -> None:
    cohort_dir = out_dir / cohort.identifier
    assemblies_dir = cohort_dir / "assemblies"
    metadata_path = cohort_dir / "metadata.tsv"

    cohort_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "ncbi-genome-download",
        cohort.group,
        "--section",
        section,
        "--formats",
        formats,
        "--species",
        cohort.species,
        "--output-folder",
        str(assemblies_dir),
        "--metadata-table",
        str(metadata_path),
        "--parallel",
        "4",
    ]
    limit = limit_override if limit_override is not None else cohort.default_limit
    if limit > 0:
        cmd.extend(["--limit", str(limit)])

    print(f"[INFO] Downloading {cohort.species} ({limit} genomes)…")
    subprocess.run(cmd, check=True)
    print(f"[INFO] Finished downloading {cohort.identifier}")


def write_manifest_and_tree(cohort: Cohort, out_dir: Path) -> None:
    cohort_dir = out_dir / cohort.identifier
    assemblies_dir = cohort_dir / "assemblies"
    manifest_path = cohort_dir / f"{cohort.identifier}.txt"
    tree_path = cohort_dir / f"{cohort.identifier}.nwk"

    fasta_paths = sorted(assemblies_dir.glob("**/*.fna.gz"))
    if not fasta_paths:
        print(
            f"[WARN] No FASTA files found for {cohort.identifier}; "
            "manifest will be empty.",
            file=sys.stderr,
        )

    manifest_path.write_text(
        "\n".join(str(path.resolve()) for path in fasta_paths) + ("\n" if fasta_paths else "")
    )
    if not tree_path.exists():
        tree_path.write_text(
            "# TODO: Replace this stub with a Newick tree covering the manifest samples.\n"
        )
    print(f"[INFO] Wrote manifest: {manifest_path}")
    print(f"[INFO] Tree placeholder: {tree_path}")


def main(argv: Optional[Iterable[str]] = None) -> None:
    args = parse_args(argv)
    ensure_tool_present()
    out_dir = args.output.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    selected = (
        [c for c in COHORTS if c.identifier in args.cohort]
        if args.cohort
        else COHORTS
    )

    for cohort in selected:
        run_download(
            cohort=cohort,
            out_dir=out_dir,
            formats=args.formats,
            section=args.section,
            limit_override=args.limit,
        )
        write_manifest_and_tree(cohort=cohort, out_dir=out_dir)

    print("[INFO] All cohorts processed.")


if __name__ == "__main__":
    main()
