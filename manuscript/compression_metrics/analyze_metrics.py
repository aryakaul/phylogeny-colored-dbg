#!/usr/bin/env python3
"""Aggregate pcDBG compression metrics for manuscript reporting.

Given one or more TSV files produced by `compute_compression_metrics`, this
script extracts summary statistics and writes a consolidated table that can be
referenced in the manuscript.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

METRIC_KEYS = (
    "cuttlefish_unique_colors",
    "pcdbg_unique_colors",
    "absolute_difference",
    "pcdbg_to_cuttlefish_ratio",
)

FILENAME_PATTERN = re.compile(
    r"(?P<batch>.+)_compression_(?P<label>.+)_k(?P<kmer>\d+)\.tsv$"
)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarise pcDBG compression metric TSV files."
    )
    parser.add_argument(
        "--metrics",
        type=Path,
        action="append",
        required=True,
        help="Path to a compression metrics TSV (may be specified multiple times).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional destination TSV for the aggregated summary (stdout if omitted).",
    )
    return parser.parse_args(argv)


def _metadata_from_path(path: Path) -> Dict[str, str]:
    match = FILENAME_PATTERN.search(path.name)
    if match:
        return {
            "batch": match.group("batch"),
            "parsimony_label": match.group("label"),
            "kmer": match.group("kmer"),
        }
    return {
        "batch": path.stem,
        "parsimony_label": "unknown",
        "kmer": "unknown",
    }


def _read_metrics(path: Path) -> Dict[str, str]:
    metrics: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            metric = row.get("metric")
            value = row.get("value")
            if not metric or value is None:
                continue
            metrics[metric] = value
    return metrics


def summarise_metrics(paths: Iterable[Path]) -> List[Dict[str, str]]:
    summaries: List[Dict[str, str]] = []
    for path in paths:
        data = _read_metrics(path)
        missing = [metric for metric in METRIC_KEYS if metric not in data]
        if missing:
            raise ValueError(
                f"{path} is missing expected metrics: {', '.join(missing)}"
            )

        summary = _metadata_from_path(path)
        for key in METRIC_KEYS:
            summary[key] = data[key]
        summaries.append(summary)
    return summaries


def write_summary(
    summaries: Sequence[Dict[str, str]], destination: Optional[Path]
) -> None:
    fieldnames = [
        "batch",
        "parsimony_label",
        "kmer",
        *METRIC_KEYS,
    ]
    if destination:
        destination.parent.mkdir(parents=True, exist_ok=True)
        handle = destination.open("w", encoding="utf-8", newline="")
    else:
        handle = sys.stdout

    try:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in summaries:
            writer.writerow(row)
    finally:
        if destination:
            handle.close()


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    summaries = summarise_metrics(args.metrics)
    write_summary(summaries, args.output)


if __name__ == "__main__":
    main()
