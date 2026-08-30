#!/usr/bin/env python3
"""
Map pcDBG boundary severity results onto a reference genome.

Produces a multi-track genome-wide visualization showing where
evolutionary boundaries cluster, annotated with genomic features.

Usage:
    python map_severity_to_reference.py \
        --severity boundary_severity.tsv \
        --db pcdbg.sqldb \
        --reference USA300.gbk \
        --output output_prefix \
        --threads 16
"""

import argparse
import csv
import logging
import os
import sqlite3
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import seaborn as sns

try:
    from Bio import SeqIO
except ImportError:
    print("BioPython required: pip install biopython")
    sys.exit(1)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Style ─────────────────────────────────────────────────────────────────────

sns.set_theme(style="white", font_scale=1.0)
plt.rcParams.update({
    "svg.fonttype": "none",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})


# ── Step 1: Extract unitig sequences ─────────────────────────────────────────

def extract_unitig_sequences(db_path, output_fasta):
    """Dump all unitig sequences from SQLite to FASTA."""
    log.info("Extracting unitig sequences from SQLite...")
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    n = 0
    with open(output_fasta, "w") as f:
        for uid, seq in cur.execute("SELECT unitig_id, sequence FROM unitigs"):
            if seq:
                f.write(f">{uid}\n{seq}\n")
                n += 1
    conn.close()
    log.info(f"  Wrote {n:,} unitig sequences to {output_fasta}")
    return n


# ── Step 2: BLAST unitigs against reference ──────────────────────────────────

def extract_reference_fasta(gbk_path, output_fasta):
    """Extract sequence from GenBank file to FASTA."""
    record = SeqIO.read(gbk_path, "genbank")
    with open(output_fasta, "w") as f:
        f.write(f">{record.id}\n{str(record.seq)}\n")
    log.info(f"  Reference: {record.id}, {len(record.seq):,} bp")
    return record


def run_blast(query_fasta, ref_fasta, output_tsv, threads=16, evalue=1e-10):
    """Build BLAST db and run blastn."""
    log.info("Building BLAST database...")
    subprocess.run(
        ["makeblastdb", "-in", ref_fasta, "-dbtype", "nucl", "-out", ref_fasta + ".db"],
        check=True, capture_output=True,
    )

    log.info(f"Running BLAST with {threads} threads...")
    cmd = [
        "blastn",
        "-query", query_fasta,
        "-db", ref_fasta + ".db",
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
        "-evalue", str(evalue),
        "-max_target_seqs", "1",
        "-dust", "no",
        "-num_threads", str(threads),
        "-out", output_tsv,
    ]
    subprocess.run(cmd, check=True)
    log.info(f"  BLAST output written to {output_tsv}")


def parse_blast_hits(blast_tsv, min_pident=90.0, min_coverage=0.5):
    """Parse BLAST hits, return dict of unitig_id → (ref_start, ref_end, midpoint)."""
    hits = {}
    n_total = 0
    n_pass = 0

    with open(blast_tsv) as f:
        for line in f:
            n_total += 1
            parts = line.strip().split("\t")
            qseqid = parts[0]
            pident = float(parts[2])
            aln_len = int(parts[3])
            sstart = int(parts[8])
            send = int(parts[9])
            qlen = int(parts[12])

            coverage = aln_len / qlen if qlen > 0 else 0

            if pident >= min_pident and coverage >= min_coverage:
                ref_start = min(sstart, send)
                ref_end = max(sstart, send)
                midpoint = (ref_start + ref_end) / 2
                # Keep best hit per unitig (first encountered with max_target_seqs=1)
                if qseqid not in hits:
                    hits[qseqid] = (ref_start, ref_end, midpoint)
                    n_pass += 1

    log.info(f"  BLAST hits: {n_total:,} total, {n_pass:,} passing filters "
             f"({100*n_pass/max(n_total,1):.1f}%)")
    return hits


# ── Step 3: Load severity data and map edges ─────────────────────────────────

def load_severity_and_map_edges(severity_tsv, unitig_hits):
    """Load severity TSV, map edges to reference positions."""
    edges_mapped = []
    edges_partial = 0
    edges_unmapped = 0
    edges_total = 0

    # Track which unitigs appear and whether they mapped, for unmapped fraction
    unitig_in_mapped_edge = set()
    unitig_in_unmapped_edge = set()

    with open(severity_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            edges_total += 1
            ua = row["unitig_a"]
            ub = row["unitig_b"]
            concordant = row["concordant"] == "True"
            try:
                tree_dist = float(row["tree_distance"])
            except (ValueError, TypeError):
                tree_dist = 0.0
            try:
                norm_dist = float(row["normalized_distance"])
            except (ValueError, TypeError):
                norm_dist = 0.0
            try:
                pvalue = float(row["pvalue"])
            except (ValueError, TypeError):
                pvalue = np.nan

            a_mapped = ua in unitig_hits
            b_mapped = ub in unitig_hits

            if a_mapped and b_mapped:
                mid_a = unitig_hits[ua][2]
                mid_b = unitig_hits[ub][2]
                edge_pos = (mid_a + mid_b) / 2
                edges_mapped.append({
                    "position": edge_pos,
                    "concordant": concordant,
                    "tree_distance": tree_dist,
                    "normalized_distance": norm_dist,
                    "pvalue": pvalue,
                })
                unitig_in_mapped_edge.update([ua, ub])
            elif a_mapped or b_mapped:
                edges_partial += 1
                # Record the mapped unitig's position for unmapped fraction tracking
                mapped_id = ua if a_mapped else ub
                unitig_in_unmapped_edge.add(mapped_id)
            else:
                edges_unmapped += 1

    log.info(f"  Edges: {edges_total:,} total, {len(edges_mapped):,} both mapped, "
             f"{edges_partial:,} partial, {edges_unmapped:,} unmapped")

    df = pd.DataFrame(edges_mapped)

    return df, unitig_hits, unitig_in_unmapped_edge


# ── Step 4: Sliding window statistics ─────────────────────────────────────────

def compute_window_stats(edges_df, unitig_hits, unmapped_edge_unitigs,
                         genome_length, window_size, step_size):
    """Compute per-window concordance/discordance statistics."""
    log.info(f"Computing window stats (window={window_size}, step={step_size})...")

    # Precompute mapped unitig positions for unmapped fraction
    mapped_unitig_positions = {}
    for uid, (rstart, rend, mid) in unitig_hits.items():
        mapped_unitig_positions[uid] = mid

    # Unmapped-edge unitig positions
    unmapped_positions = []
    for uid in unmapped_edge_unitigs:
        if uid in mapped_unitig_positions:
            unmapped_positions.append(mapped_unitig_positions[uid])
    unmapped_positions = np.array(unmapped_positions) if unmapped_positions else np.array([])

    # Edge positions as numpy arrays for fast windowing
    if len(edges_df) == 0:
        log.warning("No mapped edges to compute window stats")
        return pd.DataFrame()

    positions = edges_df["position"].values
    concordant = edges_df["concordant"].values
    tree_dists = edges_df["tree_distance"].values

    windows = []
    starts = np.arange(0, genome_length - window_size + 1, step_size)

    for w_start in starts:
        w_end = w_start + window_size

        # Edges in window
        mask = (positions >= w_start) & (positions < w_end)
        n_edges = mask.sum()

        if n_edges == 0:
            windows.append({
                "window_start": int(w_start),
                "window_end": int(w_end),
                "n_edges": 0,
                "n_concordant": 0,
                "n_discordant": 0,
                "concordance_rate": np.nan,
                "mean_tree_distance": np.nan,
                "n_unmapped_nearby": 0,
                "unmapped_fraction": 0.0,
            })
            continue

        conc = concordant[mask]
        n_conc = conc.sum()
        n_disc = n_edges - n_conc

        disc_dists = tree_dists[mask & ~concordant]
        mean_td = disc_dists.mean() if len(disc_dists) > 0 else np.nan

        # Unmapped edge unitigs nearby
        if len(unmapped_positions) > 0:
            n_unmapped = ((unmapped_positions >= w_start) & (unmapped_positions < w_end)).sum()
        else:
            n_unmapped = 0

        total_activity = n_edges + n_unmapped
        unmapped_frac = n_unmapped / total_activity if total_activity > 0 else 0.0

        windows.append({
            "window_start": int(w_start),
            "window_end": int(w_end),
            "n_edges": int(n_edges),
            "n_concordant": int(n_conc),
            "n_discordant": int(n_disc),
            "concordance_rate": n_conc / n_edges,
            "mean_tree_distance": mean_td,
            "n_unmapped_nearby": int(n_unmapped),
            "unmapped_fraction": unmapped_frac,
        })

    df = pd.DataFrame(windows)
    log.info(f"  {len(df):,} windows computed")
    return df


# ── Step 5: Parse GenBank features ────────────────────────────────────────────

MGE_KEYWORDS_INTEGRASE = ["integrase", "transposase", "recombinase", "resolvase"]
MGE_KEYWORDS_ELEMENT = [
    "phage", "prophage", "mec", "sccmec", "acme", "arca", "arcd",
    "sapi", "pathogenicity island", "insertion sequence",
    "is element", "tn916", "tn552", "tn554",
]
LANDMARK_GENES = ["orfx", "rlmn", "meca", "mecc", "arca", "arcd", "pvl", "luk"]


def classify_feature(gene_name, product):
    """Classify a CDS feature into MGE_integrase, MGE_element, landmark, or core."""
    text = f"{gene_name} {product}".lower()

    for kw in MGE_KEYWORDS_INTEGRASE:
        if kw in text:
            return "MGE_integrase"
    for kw in MGE_KEYWORDS_ELEMENT:
        if kw in text:
            return "MGE_element"
    for kw in LANDMARK_GENES:
        if kw in text:
            return "landmark"
    return "core"


def parse_genbank_features(gbk_path):
    """Extract CDS/gene features with classifications."""
    record = SeqIO.read(gbk_path, "genbank")
    features = []

    for feat in record.features:
        if feat.type not in ("CDS", "gene"):
            continue

        gene = feat.qualifiers.get("gene", [""])[0]
        product = feat.qualifiers.get("product", [""])[0]
        note = feat.qualifiers.get("note", [""])[0]

        start = int(feat.location.start)
        end = int(feat.location.end)
        category = classify_feature(gene, f"{product} {note}")

        features.append({
            "start": start,
            "end": end,
            "gene": gene,
            "product": product,
            "category": category,
        })

    df = pd.DataFrame(features)
    n_mge = len(df[df["category"].isin(["MGE_integrase", "MGE_element", "landmark"])])
    log.info(f"  Parsed {len(df):,} features ({n_mge:,} MGE-related)")
    return df


# ── Step 6: Feature proximity enrichment ──────────────────────────────────────

def compute_feature_enrichment(edges_df, features_df, genome_length, flank=2500):
    """Compute concordance rates near each non-core feature."""
    if len(edges_df) == 0 or len(features_df) == 0:
        return pd.DataFrame()

    positions = edges_df["position"].values
    concordant = edges_df["concordant"].values
    tree_dists = edges_df["tree_distance"].values

    # Background rates
    bg_concordance = concordant.mean()
    disc_mask = ~concordant
    bg_tree_dist = tree_dists[disc_mask].mean() if disc_mask.any() else np.nan

    mge_features = features_df[features_df["category"] != "core"].copy()
    results = []

    for _, feat in mge_features.iterrows():
        w_start = max(0, feat["start"] - flank)
        w_end = min(genome_length, feat["end"] + flank)

        mask = (positions >= w_start) & (positions < w_end)
        n_edges = mask.sum()

        if n_edges < 5:
            continue

        conc_rate = concordant[mask].mean()
        disc_local = tree_dists[mask & ~concordant]
        mean_td = disc_local.mean() if len(disc_local) > 0 else np.nan

        discordance_rate = 1 - conc_rate
        bg_discordance = 1 - bg_concordance
        fold_change = discordance_rate / bg_discordance if bg_discordance > 0 else np.nan

        results.append({
            "gene": feat["gene"],
            "product": feat["product"],
            "category": feat["category"],
            "start": feat["start"],
            "end": feat["end"],
            "n_edges_nearby": n_edges,
            "concordance_rate": conc_rate,
            "bg_concordance_rate": bg_concordance,
            "mean_tree_dist": mean_td,
            "bg_mean_tree_dist": bg_tree_dist,
            "fold_change_discordance": fold_change,
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        df = df.sort_values("fold_change_discordance", ascending=False)
    log.info(f"  Feature enrichment computed for {len(df):,} MGE-related features")
    return df


# ── Step 7: Visualization ─────────────────────────────────────────────────────

def plot_genome_heatmap(window_df, features_df, genome_length, output_path,
                        reference_name="Reference"):
    """Multi-track genome-wide visualization."""
    log.info("Generating genome-wide plot...")

    if len(window_df) == 0:
        log.warning("No window data to plot")
        return

    fig, axes = plt.subplots(5, 1, figsize=(16, 10), sharex=True,
                              gridspec_kw={"height_ratios": [1.5, 1.5, 1, 1, 0.8]},
                              constrained_layout=True)

    x_kb = window_df["window_start"].values / 1000
    genome_kb = genome_length / 1000

    # ── Track 1: Concordance rate heatmap ─────────────────────────────────
    ax = axes[0]
    conc = window_df["concordance_rate"].values
    bg_conc = np.nanmean(conc)

    # Create a 2D array for imshow (single row)
    conc_2d = conc.reshape(1, -1)
    vmin = max(0, bg_conc - 0.3)
    vmax = min(1, bg_conc + 0.3)

    # Use pcolormesh for correct x-axis scaling
    x_edges = np.append(x_kb, x_kb[-1] + (x_kb[1] - x_kb[0]) if len(x_kb) > 1 else x_kb[-1] + 1)
    y_edges = [0, 1]
    pcm = ax.pcolormesh(x_edges, y_edges, conc_2d,
                        cmap="RdBu", vmin=vmin, vmax=vmax, shading="flat")
    ax.set_yticks([])
    ax.set_ylabel("Concordance\nrate", fontsize=9)
    ax.set_title(f"Boundary severity mapped to {reference_name}", fontweight="bold",
                 fontsize=12, loc="left")
    cb = fig.colorbar(pcm, ax=ax, orientation="vertical", shrink=0.8, pad=0.02)
    cb.set_label("Concordance rate", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    # ── Track 2: Mean tree distance of discordant edges ──────────────────
    ax = axes[1]
    mean_td = window_df["mean_tree_distance"].values
    has_disc = ~np.isnan(mean_td)

    td_2d = mean_td.copy().reshape(1, -1)
    td_2d[0, ~has_disc] = 0  # mask NaN for display

    # Only show color where discordant edges exist
    masked = np.ma.masked_where(np.isnan(mean_td.reshape(1, -1)), td_2d)
    pcm2 = ax.pcolormesh(x_edges, y_edges, masked,
                         cmap="YlOrRd", vmin=0,
                         vmax=np.nanpercentile(mean_td[has_disc], 95) if has_disc.any() else 1,
                         shading="flat")
    ax.set_yticks([])
    ax.set_ylabel("Tree dist\n(discordant)", fontsize=9)
    cb2 = fig.colorbar(pcm2, ax=ax, orientation="vertical", shrink=0.8, pad=0.02)
    cb2.set_label("Mean tree distance", fontsize=8)
    cb2.ax.tick_params(labelsize=7)

    # ── Track 3: Edge density ─────────────────────────────────────────────
    ax = axes[2]
    ax.fill_between(x_kb, window_df["n_edges"].values, color="#4a90d9", alpha=0.6,
                    linewidth=0)
    ax.set_ylabel("Edge\ndensity", fontsize=9)
    ax.set_ylim(bottom=0)
    ax.tick_params(axis="y", labelsize=7)
    sns.despine(ax=ax)

    # ── Track 4: Unmapped fraction ────────────────────────────────────────
    ax = axes[3]
    ax.fill_between(x_kb, window_df["unmapped_fraction"].values,
                    color="#e74c3c", alpha=0.5, linewidth=0)
    ax.set_ylabel("Unmapped\nfraction", fontsize=9)
    ax.set_ylim(0, min(1, window_df["unmapped_fraction"].max() * 1.2 + 0.01))
    ax.tick_params(axis="y", labelsize=7)
    sns.despine(ax=ax)

    # ── Track 5: Gene annotations ─────────────────────────────────────────
    ax = axes[4]
    category_colors = {
        "MGE_integrase": "#e74c3c",
        "MGE_element": "#f39c12",
        "landmark": "#9b59b6",
        "core": "#cccccc",
    }

    # Only plot non-core features to reduce clutter
    mge_feats = features_df[features_df["category"] != "core"]

    for _, feat in mge_feats.iterrows():
        start_kb = feat["start"] / 1000
        width_kb = (feat["end"] - feat["start"]) / 1000
        color = category_colors.get(feat["category"], "#cccccc")
        ax.barh(0, width_kb, left=start_kb, height=0.6, color=color,
                edgecolor="none", alpha=0.8)

        # Label landmark genes
        gene = feat["gene"]
        if gene and any(kw in gene.lower() for kw in LANDMARK_GENES):
            mid_kb = start_kb + width_kb / 2
            ax.text(mid_kb, 0.8, gene, fontsize=6, ha="center", va="bottom",
                    rotation=45, style="italic")

    ax.set_ylim(-0.5, 1.5)
    ax.set_yticks([])
    ax.set_ylabel("Features", fontsize=9)
    ax.set_xlabel("Genome position (kb)", fontsize=10)
    ax.set_xlim(0, genome_kb)
    sns.despine(ax=ax, left=True)

    # Legend for feature categories
    from matplotlib.patches import Patch
    legend_patches = [
        Patch(color="#e74c3c", label="Integrase/transposase"),
        Patch(color="#f39c12", label="MGE element"),
        Patch(color="#9b59b6", label="Landmark gene"),
    ]
    ax.legend(handles=legend_patches, loc="upper right", fontsize=7,
              frameon=False, ncol=3)

    fig.savefig(output_path, format="svg", dpi=150)
    log.info(f"  Plot saved to {output_path}")
    plt.close()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Map pcDBG boundary severity to a reference genome"
    )
    parser.add_argument("--severity", required=True, help="Boundary severity TSV")
    parser.add_argument("--db", required=True, help="pcDBG SQLite database")
    parser.add_argument("--reference", required=True, help="Reference GenBank file")
    parser.add_argument("--output", required=True, help="Output prefix")
    parser.add_argument("--window", type=int, default=5000, help="Window size (bp)")
    parser.add_argument("--step", type=int, default=1000, help="Step size (bp)")
    parser.add_argument("--threads", type=int, default=16, help="BLAST threads")
    parser.add_argument("--evalue", type=float, default=1e-10, help="BLAST evalue")
    parser.add_argument("--min-pident", type=float, default=90.0,
                        help="Min percent identity for BLAST hits")
    parser.add_argument("--min-coverage", type=float, default=0.5,
                        help="Min query coverage for BLAST hits")
    parser.add_argument("--flank", type=int, default=2500,
                        help="Flank size for feature enrichment (bp)")
    args = parser.parse_args()

    output_dir = os.path.dirname(args.output) or "."
    os.makedirs(output_dir, exist_ok=True)

    # ── Step 1: Extract unitig sequences ──────────────────────────────────
    unitig_fasta = args.output + "_unitigs.fasta"
    if not os.path.exists(unitig_fasta):
        extract_unitig_sequences(args.db, unitig_fasta)
    else:
        log.info(f"Using existing unitig FASTA: {unitig_fasta}")

    # ── Step 2: BLAST ─────────────────────────────────────────────────────
    log.info("Parsing reference GenBank...")
    ref_fasta = args.output + "_reference.fasta"
    record = extract_reference_fasta(args.reference, ref_fasta)
    genome_length = len(record.seq)
    reference_name = record.id
    if record.description and record.description != record.id:
        reference_name = record.description[:80]

    blast_output = args.output + "_blast_hits.tsv"
    if not os.path.exists(blast_output):
        run_blast(unitig_fasta, ref_fasta, blast_output,
                  threads=args.threads, evalue=args.evalue)
    else:
        log.info(f"Using existing BLAST output: {blast_output}")

    log.info("Parsing BLAST hits...")
    unitig_hits = parse_blast_hits(blast_output,
                                   min_pident=args.min_pident,
                                   min_coverage=args.min_coverage)
    log.info(f"  {len(unitig_hits):,} unitigs mapped to reference")

    # ── Step 3: Map edges to reference ────────────────────────────────────
    log.info("Loading severity data and mapping edges...")
    edges_df, _, unmapped_unitigs = load_severity_and_map_edges(
        args.severity, unitig_hits
    )

    # Save unitig map
    unitig_map_path = args.output + "_unitig_map.tsv"
    with open(unitig_map_path, "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["unitig_id", "ref_start", "ref_end", "ref_midpoint"])
        for uid, (rs, re, mid) in sorted(unitig_hits.items(),
                                          key=lambda x: x[1][2]):
            w.writerow([uid, rs, re, f"{mid:.0f}"])
    log.info(f"  Unitig map written to {unitig_map_path}")

    # ── Step 4: Window stats ──────────────────────────────────────────────
    window_df = compute_window_stats(
        edges_df, unitig_hits, unmapped_unitigs,
        genome_length, args.window, args.step,
    )
    window_path = args.output + "_window_stats.tsv"
    window_df.to_csv(window_path, sep="\t", index=False, float_format="%.6f")
    log.info(f"  Window stats written to {window_path}")

    # Print summary
    if len(window_df) > 0:
        mean_conc = window_df["concordance_rate"].mean()
        log.info(f"  Genome-wide mean concordance: {mean_conc:.4f}")
        high_disc = window_df[window_df["concordance_rate"] < mean_conc - 0.15]
        log.info(f"  Windows with concordance >15% below mean: {len(high_disc):,}")

    # ── Step 5: Parse features ────────────────────────────────────────────
    log.info("Parsing GenBank features...")
    features_df = parse_genbank_features(args.reference)

    # ── Step 6: Feature enrichment ────────────────────────────────────────
    log.info("Computing feature proximity enrichment...")
    enrichment_df = compute_feature_enrichment(
        edges_df, features_df, genome_length, flank=args.flank,
    )
    enrichment_path = args.output + "_feature_enrichment.tsv"
    if len(enrichment_df) > 0:
        enrichment_df.to_csv(enrichment_path, sep="\t", index=False,
                             float_format="%.4f")
        log.info(f"  Feature enrichment written to {enrichment_path}")

        # Print top features by discordance enrichment
        log.info("  Top features by discordance enrichment:")
        for _, row in enrichment_df.head(20).iterrows():
            log.info(f"    {row['gene']:15s} {row['product'][:40]:40s} "
                     f"fold={row['fold_change_discordance']:.2f} "
                     f"conc={row['concordance_rate']:.3f} "
                     f"n={row['n_edges_nearby']}")

    # ── Step 7: Plot ──────────────────────────────────────────────────────
    plot_path = args.output + "_genome_heatmap.svg"
    plot_genome_heatmap(window_df, features_df, genome_length, plot_path,
                        reference_name=reference_name)

    log.info("Done.")


if __name__ == "__main__":
    main()
