#!/usr/bin/env python3
"""Download high-quality genome cohorts for pcDBG experiments.

This script ensures that we retrieve at least 500 high-quality genomes for
each target genus. It first performs a dry-run metadata sweep with
`ncbi-genome-download`, ranks assemblies by curation status and assembly
level, then downloads the selected cohort and writes manifests compatible
with the pcDBG Snakemake workflow.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple
from urllib.error import URLError
from urllib.request import urlopen


HIGH_QUALITY_ASSEMBLY_LEVELS = ("complete", "chromosome")
ASSEMBLY_LEVEL_ARG = ",".join(HIGH_QUALITY_ASSEMBLY_LEVELS)
HIGH_QUALITY_CATEGORIES = ("reference", "representative", "na")
REFSEQ_CATEGORY_ARG = ",".join(HIGH_QUALITY_CATEGORIES)

ASSEMBLY_SUMMARY_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes"
ASSEMBLY_LEVEL_CANONICAL = {
    "complete": "complete_genome",
    "complete_genome": "complete_genome",
    "chromosome": "chromosome",
    "scaffold": "scaffold",
    "contig": "contig",
}

CATEGORY_CANONICAL = {
    "reference": "reference_genome",
    "reference_genome": "reference_genome",
    "representative": "representative_genome",
    "representative_genome": "representative_genome",
    "na": "na",
}

ALLOWED_ASSEMBLY_LEVELS = {"complete_genome", "chromosome"}
CATEGORY_PRIORITY = {"reference_genome": 0, "representative_genome": 1, "na": 2, "": 3}
ASSEMBLY_PRIORITY = {"complete_genome": 0, "chromosome": 1, "scaffold": 2, "contig": 3, "": 4}


@dataclass(frozen=True)
class Cohort:
    identifier: str
    group: str
    genus: str
    target_count: int = 500


COHORTS: List[Cohort] = [
    Cohort(identifier="listeria", group="bacteria", genus="Listeria"),
    Cohort(identifier="mycobacterium", group="bacteria", genus="Mycobacterium"),
    Cohort(identifier="escherichia", group="bacteria", genus="Escherichia"),
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
        choices=("refseq", "genbank"),
        help="NCBI section to download from (default: refseq).",
    )
    parser.add_argument(
        "--min-genomes",
        type=int,
        help="Override the minimum number of genomes to fetch per cohort (default: 500).",
    )
    parser.add_argument(
        "--cohort",
        action="append",
        choices=[c.identifier for c in COHORTS],
        help="Download only the specified cohort(s). Can be used multiple times.",
    )
    parser.add_argument(
        "--parallel",
        type=int,
        default=4,
        help="Number of parallel downloads to run (passed to ncbi-genome-download).",
    )
    return parser.parse_args(argv)


def ensure_tool_present() -> None:
    if shutil.which("ncbi-genome-download") is None:
        raise RuntimeError(
            "ncbi-genome-download is not available on PATH. "
            "Install it via `pip install ncbi-genome-download`."
        )


def ensure_fasta_requested(formats: str) -> None:
    requested = {fmt.strip().lower() for fmt in formats.split(",") if fmt.strip()}
    if "fasta" not in requested:
        raise RuntimeError(
            "The download formats must include 'fasta' so manifests can reference FASTA files."
        )


def canonicalize_assembly_level(value: Optional[str]) -> str:
    if not value:
        return ""
    cleaned = value.strip().lower().replace(" ", "_")
    return ASSEMBLY_LEVEL_CANONICAL.get(cleaned, cleaned)


def canonicalize_refseq_category(value: Optional[str]) -> str:
    if not value:
        return ""
    cleaned = value.strip().lower().replace(" ", "_")
    return CATEGORY_CANONICAL.get(cleaned, cleaned)


def is_genus_match(organism_name: str, genus: str) -> bool:
    if not organism_name:
        return False
    tokens = organism_name.strip().split()
    if not tokens:
        return False
    return tokens[0].lower() == genus.lower()


def fetch_genus_metadata(
    section: str,
    group: str,
    genus: str,
) -> Tuple[List[dict], List[str]]:
    url = f"{ASSEMBLY_SUMMARY_BASE}/{section}/{group}/assembly_summary.txt"
    header: Optional[List[str]] = None
    rows: List[dict] = []

    try:
        with urlopen(url) as response:
            for raw_line in response:
                line = raw_line.decode("utf-8").rstrip("\n")
                if not line:
                    continue
                if line.startswith("#"):
                    header = line.lstrip("#").strip().split("\t")
                    continue
                if not header:
                    continue
                parts = line.split("\t")
                if len(parts) < len(header):
                    parts.extend([""] * (len(header) - len(parts)))
                row = {header[i]: parts[i] for i in range(len(header))}
                if not is_genus_match(row.get("organism_name", ""), genus):
                    continue
                if (row.get("version_status") or "").lower() != "latest":
                    continue
                rows.append(row)
    except URLError as exc:
        raise RuntimeError(f"Failed to fetch assembly summary from {url}: {exc}") from exc

    if not rows:
        raise RuntimeError(
            f"No assemblies found for genus '{genus}' in section '{section}' group '{group}'."
        )

    return rows, header or []


def release_date_score(row: dict) -> int:
    for key in ("seq_rel_date", "release_date"):
        raw = (row.get(key) or "").strip()
        if not raw:
            continue
        try:
            parsed = dt.datetime.strptime(raw, "%Y-%m-%d").date()
            return -int(parsed.strftime("%Y%m%d"))
        except ValueError:
            continue
    return 0


def priority_key(row: dict) -> Tuple[int, int, int, str]:
    category = canonicalize_refseq_category(row.get("refseq_category"))
    assembly = canonicalize_assembly_level(row.get("assembly_level"))
    return (
        CATEGORY_PRIORITY.get(category, CATEGORY_PRIORITY[""]),
        ASSEMBLY_PRIORITY.get(assembly, ASSEMBLY_PRIORITY[""]),
        release_date_score(row),
        row.get("assembly_accession", ""),
    )


def select_high_quality(rows: Sequence[dict], minimum: int) -> List[dict]:
    candidates = [
        row
        for row in rows
        if canonicalize_assembly_level(row.get("assembly_level")) in ALLOWED_ASSEMBLY_LEVELS
    ]
    if len(candidates) < minimum:
        raise RuntimeError(
            f"Only {len(candidates)} assemblies meet the quality criteria (complete/chromosome). "
            f"Need at least {minimum} genomes."
        )

    ordered = sorted(candidates, key=priority_key)
    selected: List[dict] = []
    seen: set[str] = set()
    for row in ordered:
        accession = (row.get("assembly_accession") or "").strip()
        if not accession or accession in seen:
            continue
        selected.append(row)
        seen.add(accession)
        if len(selected) == minimum:
            break

    if len(selected) < minimum:
        raise RuntimeError(
            f"After removing duplicates only {len(selected)} assemblies were available (need {minimum})."
        )
    return selected


def write_metadata(fieldnames: Sequence[str], rows: Sequence[dict], destination: Path) -> None:
    resolved_fields: List[str]
    if fieldnames:
        resolved_fields = list(fieldnames)
    elif rows:
        resolved_fields = sorted(rows[0].keys())
    else:
        raise RuntimeError("Cannot write metadata without field names or rows.")

    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=resolved_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def cache_metadata_table(
    section: str,
    group: str,
    genus: str,
    destination: Path,
) -> Tuple[List[dict], List[str]]:
    rows, fieldnames = fetch_genus_metadata(section=section, group=group, genus=genus)
    write_metadata(fieldnames, rows, destination)
    return rows, fieldnames


def run_download(
    cohort: Cohort,
    out_dir: Path,
    formats: str,
    section: str,
    min_genomes_override: Optional[int],
    parallel: int,
) -> List[Path]:
    target = cohort.target_count
    if min_genomes_override is not None and min_genomes_override > target:
        target = min_genomes_override

    cohort_dir = out_dir / cohort.identifier
    assemblies_dir = cohort_dir / "assemblies"
    metadata_candidates = cohort_dir / "metadata_candidates.tsv"
    selected_accessions_path = cohort_dir / "selected_accessions.txt"
    metadata_subset_path = cohort_dir / "metadata.tsv"

    cohort_dir.mkdir(parents=True, exist_ok=True)
    if assemblies_dir.exists():
        shutil.rmtree(assemblies_dir)
    assemblies_dir.mkdir(parents=True, exist_ok=True)

    for path in (metadata_candidates, selected_accessions_path, metadata_subset_path):
        if path.exists():
            path.unlink()

    print(
        f"[INFO] Resolving catalog for genus {cohort.genus} in {cohort.group} "
        f"(targeting {target} genomes)…"
    )
    metadata_rows, fieldnames = cache_metadata_table(
        section=section,
        group=cohort.group,
        genus=cohort.genus,
        destination=metadata_candidates,
    )
    selected_rows = select_high_quality(metadata_rows, target)

    accessions = [(row.get("assembly_accession") or "").strip() for row in selected_rows]
    selected_accessions_path.write_text("\n".join(accessions) + "\n")
    write_metadata(fieldnames, selected_rows, metadata_subset_path)

    download_cmd = [
        "ncbi-genome-download",
        cohort.group,
        "--section",
        section,
        "--formats",
        formats,
        "--assembly-accessions",
        str(selected_accessions_path),
        "--output-folder",
        str(assemblies_dir),
        "--assembly-levels",
        ASSEMBLY_LEVEL_ARG,
        "--refseq-categories",
        REFSEQ_CATEGORY_ARG,
        "--parallel",
        str(parallel),
    ]

    print(f"[INFO] Downloading {len(accessions)} assemblies for genus {cohort.genus}…")
    subprocess.run(download_cmd, check=True)

    fasta_paths = sorted(assemblies_dir.rglob("*.fna.gz"))
    if not fasta_paths:
        fasta_paths = sorted(assemblies_dir.rglob("*.fna"))

    if len(fasta_paths) < target:
        raise RuntimeError(
            f"Expected at least {target} FASTA files for genus {cohort.genus}, "
            f"but found {len(fasta_paths)} under {assemblies_dir}."
        )

    print(
        f"[INFO] Finished downloading {len(fasta_paths)} FASTA files for genus {cohort.genus}."
    )
    return fasta_paths


def write_manifest_and_tree(
    cohort: Cohort,
    out_dir: Path,
    fasta_paths: Sequence[Path],
) -> None:
    cohort_dir = out_dir / cohort.identifier
    assemblies_dir = cohort_dir / "assemblies"
    manifest_path = cohort_dir / f"{cohort.identifier}.txt"
    tree_path = cohort_dir / f"{cohort.identifier}.nwk"

    if not fasta_paths:
        print(
            f"[WARN] No FASTA files discovered for {cohort.identifier}; manifest will be empty.",
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
    print(f"[INFO] Assemblies directory: {assemblies_dir}")


def main(argv: Optional[Iterable[str]] = None) -> None:
    args = parse_args(argv)
    ensure_tool_present()
    ensure_fasta_requested(args.formats)

    if args.min_genomes is not None and args.min_genomes < 1:
        raise RuntimeError("--min-genomes must be a positive integer.")
    if args.parallel < 1:
        raise RuntimeError("--parallel must be a positive integer.")

    out_dir = args.output.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    selected = (
        [c for c in COHORTS if c.identifier in args.cohort]
        if args.cohort
        else COHORTS
    )

    for cohort in selected:
        fasta_paths = run_download(
            cohort=cohort,
            out_dir=out_dir,
            formats=args.formats,
            section=args.section,
            min_genomes_override=args.min_genomes,
            parallel=args.parallel,
        )
        write_manifest_and_tree(cohort=cohort, out_dir=out_dir, fasta_paths=fasta_paths)

    print("[INFO] All cohorts processed.")


if __name__ == "__main__":
    main()
