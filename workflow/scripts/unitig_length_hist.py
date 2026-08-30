#!/usr/bin/env python3
"""
Unitig Length Histogram from GFA1 (compressed de Bruijn graph)
Usage: python unitig_length_histogram.py <input.gfa> [output.png]
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def parse_gfa1_unitig_lengths(gfa_path: str) -> list[int]:
    """
    Parse a GFA1 file and extract unitig lengths from Segment (S) lines.
    Handles both explicit LN:i: tags and inline sequence length.
    """
    lengths = []

    with open(gfa_path, "r") as fh:
        for line in fh:
            if not line.startswith("S\t"):
                continue

            parts = line.rstrip("\n").split("\t")
            # parts[0] = 'S', parts[1] = name, parts[2] = sequence or '*'
            if len(parts) < 3:
                continue

            length = None

            # Prefer explicit LN:i: tag (O(1) per line, avoids storing sequence)
            for tag in parts[3:]:
                if tag.startswith("LN:i:"):
                    length = int(tag[5:])
                    break

            # Fall back to sequence length
            if length is None:
                seq = parts[2]
                if seq != "*":
                    length = len(seq)

            if length is not None and length > 0:
                lengths.append(length)

    return lengths


def plot_histogram(lengths: list[int], output_path: str | None = None) -> None:
    lengths_arr = np.array(lengths, dtype=np.int64)
    mean_len = lengths_arr.mean()
    n = len(lengths_arr)

    # ── figure setup ──────────────────────────────────────────────────────────
    sns.set_theme(style="white")
    fig, ax = plt.subplots(figsize=(8, 4))

    # Sturges-capped bin count – sensible for highly skewed assembly data
    n_bins = min(int(np.ceil(np.log2(n))) + 1, 120)

    # Log-scale x-axis handles the typical power-law length distribution
    log_lengths = np.log10(lengths_arr)
    bin_edges = np.linspace(log_lengths.min(), log_lengths.max(), n_bins + 1)

    sns.histplot(
        log_lengths,
        bins=bin_edges,
        log_scale=True,
        kde=True,
        color="#4C72B0",
        edgecolor="white",
        linewidth=0.4,
        alpha=0.85,
        ax=ax,
    )

    # ── mean vertical line ────────────────────────────────────────────────────
    mean_log = np.log10(mean_len)
    ax.axvline(
        mean_log,
        color="#DD4444",
        linewidth=2.2,
        linestyle="--",
        label=f"Mean = {mean_len:,.1f} bp",
        zorder=5,
    )

    # ── axes labels & formatting ───────────────────────────────────────────────
    # Replace log10 ticks with human-readable bp values
    tick_vals = np.arange(
        np.floor(log_lengths.min()), np.ceil(log_lengths.max()) + 1, 1
    )
    ax.set_xticks(tick_vals)
    ax.set_xticklabels([f"{10**v:,.0f}" for v in tick_vals], rotation=30, ha="right")

    ax.set_xlabel("Unitig Length (bp, log scale)", labelpad=10)
    ax.set_ylabel("Count", labelpad=10)
    # ax.set_title(
        # f"Unitig Length Distribution  (n = {n:,})", fontsize=15, fontweight="bold", pad=14
    # )

    # Summary stats annotation
    stats_text = (
        f"Median : {np.median(lengths_arr):,.0f} bp\n"
        # f"N50    : {compute_n50(lengths_arr):,.0f} bp\n"
        f"Max    : {lengths_arr.max():,.0f} bp"
    )
    ax.text(
        0.97, 0.95, stats_text,
        transform=ax.transAxes,
        va="top", ha="right",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.8, edgecolor="#cccccc"),
        # family="monospace",
    )

    # ax.legend(fontsize=11)
    fig.tight_layout()
    sns.despine()

    if output_path:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"Saved → {output_path}")
    else:
        plt.show()


def compute_n50(lengths: np.ndarray) -> int:
    sorted_desc = np.sort(lengths)[::-1]
    cumsum = np.cumsum(sorted_desc)
    return int(sorted_desc[np.searchsorted(cumsum, cumsum[-1] / 2)])


# ── entry point ───────────────────────────────────────────────────────────────
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python unitig_length_histogram.py <input.gfa> [output.png]")
        sys.exit(1)

    gfa_file = sys.argv[1]
    out_file = sys.argv[2] if len(sys.argv) > 2 else None

    print(f"Parsing {gfa_file} …")
    lengths = parse_gfa1_unitig_lengths(gfa_file)

    if not lengths:
        print("No unitig lengths found. Check that the file contains S-lines.")
        sys.exit(1)

    print(f"Found {len(lengths):,} unitigs. Plotting …")
    plot_histogram(lengths, out_file)
