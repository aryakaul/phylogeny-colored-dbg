#!/usr/bin/env python3
"""
Map pcDBG boundary severity to a reference genome — heatmap version.

5-track genome-wide visualization:
  Track 1: Concordance rate         (diverging heatmap: blue=concordant, red=discordant)
  Track 2: Discordant edge density  (sequential heatmap: white→red, raw counts)
  Track 3: Pangenome divergence     (sequential heatmap: white→purple, unmapped+structural)
  Track 4: Boundary tree distance   (sequential heatmap: white→orange, severity)
  Track 5: MGE gene annotations     (colored bars with labels)

Usage:
    python map_severity_to_reference.py \
        --severity boundary_severity.tsv \
        --db pcdbg.sqldb \
        --reference USA300.gbk \
        --output output_prefix
"""

import argparse
import csv
import logging
import os
import sqlite3
import subprocess
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import TwoSlopeNorm
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

sns.set_theme(style="white", font_scale=1.0)
plt.rcParams.update({
    "svg.fonttype": "none",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})


# ── Unitig extraction ─────────────────────────────────────────────────────────

def extract_unitig_sequences(db_path, output_fasta):
    log.info("Extracting unitig sequences from SQLite...")
    conn = sqlite3.connect(db_path)
    n = 0
    with open(output_fasta, "w") as f:
        for uid, seq in conn.execute("SELECT unitig_id, sequence FROM unitigs"):
            if seq:
                f.write(f">{uid}\n{seq}\n")
                n += 1
    conn.close()
    log.info(f"  Wrote {n:,} unitig sequences")
    return n


# ── BLAST ─────────────────────────────────────────────────────────────────────

def extract_reference_fasta(gbk_path, output_fasta):
    record = SeqIO.read(gbk_path, "genbank")
    with open(output_fasta, "w") as f:
        f.write(f">{record.id}\n{str(record.seq)}\n")
    log.info(f"  Reference: {record.id}, {len(record.seq):,} bp")
    return record


def run_blast(query_fasta, ref_fasta, output_tsv, threads=16, evalue=1e-10):
    db_prefix = ref_fasta + ".db"
    log.info("Building BLAST database...")
    subprocess.run(
        ["makeblastdb", "-in", ref_fasta, "-dbtype", "nucl", "-out", db_prefix],
        check=True, capture_output=True,
    )
    log.info(f"Running BLAST with {threads} threads...")
    subprocess.run([
        "blastn", "-query", query_fasta, "-db", db_prefix,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                   "qstart qend sstart send evalue bitscore qlen slen",
        "-evalue", str(evalue),
        "-max_target_seqs", "1",
        "-dust", "no",
        "-num_threads", str(threads),
        "-out", output_tsv,
    ], check=True)
    log.info(f"  BLAST done → {output_tsv}")


def parse_blast_hits(blast_tsv, min_pident=90.0, min_coverage=0.5):
    hits = {}
    n_total = n_pass = 0
    with open(blast_tsv) as f:
        for line in f:
            n_total += 1
            p = line.strip().split("\t")
            qid, pident, alen = p[0], float(p[2]), int(p[3])
            sstart, send, qlen = int(p[8]), int(p[9]), int(p[12])
            if pident >= min_pident and (alen / max(qlen, 1)) >= min_coverage and qid not in hits:
                rs, re = min(sstart, send), max(sstart, send)
                hits[qid] = (rs, re, (rs + re) / 2)
                n_pass += 1
    log.info(f"  BLAST: {n_total:,} hits, {n_pass:,} passing ({100*n_pass/max(n_total,1):.1f}%)")
    return hits


# ── Edge classification ───────────────────────────────────────────────────────

def classify_edges(severity_tsv, unitig_hits, max_ref_dist=10000):
    local_concordant = []
    local_discordant = []
    structural_positions = []
    unmapped_positions = []

    n_total = n_struct = n_partial = n_both_unmapped = 0

    with open(severity_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            n_total += 1
            ua, ub = row["unitig_a"], row["unitig_b"]
            concordant = row["concordant"] == "True"
            try:
                tree_dist = float(row["tree_distance"])
            except (ValueError, TypeError):
                tree_dist = 0.0

            a_hit, b_hit = unitig_hits.get(ua), unitig_hits.get(ub)

            if a_hit and b_hit:
                mid_a, mid_b = a_hit[2], b_hit[2]
                if abs(mid_a - mid_b) <= max_ref_dist:
                    pos = (mid_a + mid_b) / 2
                    if concordant:
                        local_concordant.append(pos)
                    else:
                        local_discordant.append((pos, tree_dist))
                else:
                    n_struct += 1
                    structural_positions.extend([mid_a, mid_b])
            elif a_hit or b_hit:
                n_partial += 1
                unmapped_positions.append((a_hit or b_hit)[2])
            else:
                n_both_unmapped += 1

    log.info(f"  Edges ({n_total:,}): {len(local_concordant):,} local-conc, "
             f"{len(local_discordant):,} local-disc, {n_struct:,} structural, "
             f"{n_partial:,} partial-unmapped, {n_both_unmapped:,} fully-unmapped")

    disc_arr = (np.array(local_discordant, dtype=[("pos", np.float64), ("td", np.float64)])
                if local_discordant
                else np.array([], dtype=[("pos", np.float64), ("td", np.float64)]))

    return (
        np.array(local_concordant, dtype=np.float64),
        disc_arr,
        np.array(structural_positions, dtype=np.float64),
        np.array(unmapped_positions, dtype=np.float64),
    )


# ── Window statistics ─────────────────────────────────────────────────────────

def compute_window_stats(local_conc, local_disc, structural, unmapped,
                         genome_length, window_size, step_size):
    log.info(f"Windows: size={window_size}, step={step_size}...")

    disc_pos = local_disc["pos"] if len(local_disc) > 0 else np.array([])
    disc_td = local_disc["td"] if len(local_disc) > 0 else np.array([])

    starts = np.arange(0, genome_length - window_size + 1, step_size)
    rows = []

    for ws in starts:
        # Cap the window end at the actual genome length
        we = min(ws + window_size, genome_length)

        n_conc = int(np.sum((local_conc >= ws) & (local_conc < we)))

        disc_mask = (disc_pos >= ws) & (disc_pos < we)
        n_disc = int(np.sum(disc_mask))

        n_struct = int(np.sum((structural >= ws) & (structural < we)))
        n_unmap = int(np.sum((unmapped >= ws) & (unmapped < we)))

        n_local = n_conc + n_disc
        conc_rate = n_conc / n_local if n_local > 0 else np.nan

        # Mean and max tree distance of discordant edges
        if n_disc > 0:
            local_td = disc_td[disc_mask]
            mean_td = float(np.mean(local_td))
            max_td = float(np.max(local_td))
        else:
            mean_td = np.nan
            max_td = np.nan

        divergence = n_struct + n_unmap

        rows.append({
            "window_start": int(ws),
            "window_end": int(we),
            "n_concordant": n_conc,
            "n_discordant": n_disc,
            "n_local_edges": n_local,
            "concordance_rate": conc_rate,
            "n_divergence": divergence,
            "mean_tree_distance": mean_td,
            "max_tree_distance": max_td,
        })

    df = pd.DataFrame(rows)
    log.info(f"  {len(df):,} windows computed")
    return df


# ── GenBank features ──────────────────────────────────────────────────────────

MGE_KW_INTEGRASE = ["integrase", "transposase", "recombinase", "resolvase"]
MGE_KW_ELEMENT = [
    "phage", "prophage", "mec", "sccmec", "acme", "arca", "arcd",
    "sapi", "pathogenicity island", "insertion sequence",
    "is element", "is256", "is431", "is1272",
    "tn916", "tn552", "tn554",
]
LANDMARK_GENES = ["orfx", "rlmn", "meca", "mecc", "arca", "arcd",
                  "pvl", "luks", "lukf", "lukd", "luke", "luk"]


def classify_feature(gene, product):
    text = f"{gene} {product}".lower()
    for kw in MGE_KW_INTEGRASE:
        if kw in text:
            return "MGE_integrase"
    for kw in MGE_KW_ELEMENT:
        if kw in text:
            return "MGE_element"
    for kw in LANDMARK_GENES:
        if kw in text:
            return "landmark"
    return "core"


def parse_features(gbk_path):
    record = SeqIO.read(gbk_path, "genbank")
    feats = []
    for ft in record.features:
        if ft.type not in ("CDS", "gene"):
            continue
        gene = ft.qualifiers.get("gene", [""])[0]
        product = ft.qualifiers.get("product", [""])[0]
        note = ft.qualifiers.get("note", [""])[0]
        cat = classify_feature(gene, f"{product} {note}")
        feats.append({
            "start": int(ft.location.start),
            "end": int(ft.location.end),
            "gene": gene,
            "product": product,
            "category": cat,
        })
    df = pd.DataFrame(feats)
    n_mge = len(df[df["category"] != "core"])
    log.info(f"  {len(df):,} features ({n_mge:,} MGE-related)")
    return df


# ── Feature enrichment ────────────────────────────────────────────────────────

def compute_enrichment(window_df, features_df, genome_length, step_size, flank=2500):
    if len(window_df) == 0 or len(features_df) == 0:
        return pd.DataFrame()

    total_disc = window_df["n_discordant"].sum()
    bg_density = total_disc / genome_length

    mge = features_df[features_df["category"] != "core"].copy()
    results = []

    for _, feat in mge.iterrows():
        ws = max(0, feat["start"] - flank)
        we = min(genome_length, feat["end"] + flank)
        region_len = we - ws

        mask = (window_df["window_start"] >= ws - step_size) & (window_df["window_start"] < we)
        if mask.sum() == 0:
            continue

        local_disc = window_df.loc[mask, "n_discordant"].sum()
        local_density = local_disc / region_len if region_len > 0 else 0
        fold = local_density / bg_density if bg_density > 0 else np.nan

        n_div = window_df.loc[mask, "n_divergence"].sum()

        results.append({
            "gene": feat["gene"],
            "product": feat["product"],
            "category": feat["category"],
            "start": feat["start"],
            "end": feat["end"],
            "discordant_edges": int(local_disc),
            "fold_enrichment": fold,
            "divergence_edges": int(n_div),
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        df = df.sort_values("fold_enrichment", ascending=False)
    log.info(f"  Enrichment for {len(df):,} MGE features")
    return df


# ── Visualization ─────────────────────────────────────────────────────────────

def _smooth(arr, k=5):
    """Simple moving average, preserving length."""
    if len(arr) < k:
        return arr
    kernel = np.ones(k) / k
    return np.convolve(arr, kernel, mode="same")


def _make_heatmap_strip(ax, x_edges, values, cmap, vmin, vmax, label,
                        norm=None, nan_color="#f0f0f0"):
    """Draw a single-row pcolormesh heatmap strip."""
    data = np.ma.masked_invalid(values.reshape(1, -1))
    kwargs = dict(cmap=cmap, shading="flat", rasterized=True)
    if norm is not None:
        kwargs["norm"] = norm
    else:
        kwargs["vmin"] = vmin
        kwargs["vmax"] = vmax
    pcm = ax.pcolormesh(x_edges, [0, 1], data, **kwargs)
    ax.set_facecolor(nan_color)
    ax.set_yticks([])
    ax.set_ylabel(label, fontsize=9, rotation=0, ha="right", va="center",
                  labelpad=10)
    return pcm


def plot_genome(window_df, features_df, genome_length, output_path,
                reference_name="Reference"):
    log.info("Generating genome heatmap...")

    if len(window_df) == 0:
        log.warning("No data to plot")
        return

    fig, axes = plt.subplots(
        5, 1, figsize=(18, 7.5), sharex=True,
        gridspec_kw={"height_ratios": [1, 1, 1, 1, 0.6],
                     "hspace": 0.08},
    )
    fig.subplots_adjust(right=0.88)

    x_kb = window_df["window_start"].values / 1000
    dx = x_kb[1] - x_kb[0] if len(x_kb) > 1 else 1
    genome_kb = genome_length / 1000
    
    # Calculate edges, then force the last edge to snap exactly to genome_kb
    x_edges = np.append(x_kb, x_kb[-1] + dx)
    x_edges[-1] = genome_kb

    # ── Track 1: Concordance rate (diverging heatmap) ─────────────────────
    # Blue = co-inherited block. Red = evolutionary boundary.
    # Centered on genome-wide mean so colors show deviation from background.
    ax = axes[0]
    conc_rate = window_df["concordance_rate"].values
    bg_conc = np.nanmean(conc_rate)

    # Smooth for visual clarity
    conc_smooth = _smooth(conc_rate, k=3)

    # TwoSlopeNorm centers the colormap at the background rate
    vmin_c = max(0, bg_conc - 0.4)
    vmax_c = min(1, bg_conc + 0.4)
    norm_c = TwoSlopeNorm(vmin=vmin_c, vcenter=bg_conc, vmax=vmax_c)

    pcm1 = _make_heatmap_strip(ax, x_edges, conc_smooth, "RdBu", vmin_c, vmax_c,
                                "Concordance\nrate", norm=norm_c)
    cb1 = fig.colorbar(pcm1, ax=ax, orientation="vertical", shrink=0.8,
                       pad=0.015, aspect=8)
    cb1.set_label(f"Rate (bg={bg_conc:.2f})", fontsize=7)
    cb1.ax.tick_params(labelsize=6)

    ax.set_title(f"Boundary severity mapped to {reference_name}",
                 fontweight="bold", fontsize=11, loc="left", pad=8)

    # ── Track 2: Discordant edge density (sequential heatmap) ─────────────
    # White = no boundaries. Dark red = dense evolutionary boundaries.
    # Raw counts preserve sharp spikes that rates dilute.
    ax = axes[1]
    disc = window_df["n_discordant"].values.astype(float)
    disc_smooth = _smooth(disc, k=3)

    p95 = np.percentile(disc_smooth[disc_smooth > 0], 95) if (disc_smooth > 0).any() else 1
    pcm2 = _make_heatmap_strip(ax, x_edges, disc_smooth, "Reds", 0, p95,
                                "Discordant\nedge density")
    cb2 = fig.colorbar(pcm2, ax=ax, orientation="vertical", shrink=0.8,
                       pad=0.015, aspect=8)
    cb2.set_label("Count / window", fontsize=7)
    cb2.ax.tick_params(labelsize=6)

    # ── Track 3: Pangenome divergence (sequential heatmap) ────────────────
    # White = reference-like. Dark purple = high structural/unmapped content.
    # Peaks mark where accessory content inserts relative to the reference.
    ax = axes[2]
    div = window_df["n_divergence"].values.astype(float)
    div_smooth = _smooth(div, k=3)

    p95_div = np.percentile(div_smooth[div_smooth > 0], 95) if (div_smooth > 0).any() else 1
    pcm3 = _make_heatmap_strip(ax, x_edges, div_smooth, "Purples", 0, p95_div,
                                "Pangenome\ndivergence")
    cb3 = fig.colorbar(pcm3, ax=ax, orientation="vertical", shrink=0.8,
                       pad=0.015, aspect=8)
    cb3.set_label("Unmapped + structural", fontsize=7)
    cb3.ax.tick_params(labelsize=6)

    # ── Track 4: Mean tree distance at discordant edges ───────────────────
    # White/pale = mild boundaries (nearby branches). Orange/red = severe.
    # Based on overall results, most should be mild — spatial variation is
    # the interesting signal here.
    ax = axes[3]
    mean_td = window_df["mean_tree_distance"].values

    # Smooth, preserving NaN pattern
    td_filled = np.where(np.isnan(mean_td), 0, mean_td)
    td_smooth = _smooth(td_filled, k=3)
    # Re-mask windows that originally had no discordant edges
    td_smooth[window_df["n_discordant"].values == 0] = np.nan

    p95_td = np.nanpercentile(td_smooth, 95) if np.any(~np.isnan(td_smooth)) else 1
    pcm4 = _make_heatmap_strip(ax, x_edges, td_smooth, "YlOrRd", 0, p95_td,
                                "Boundary\nseverity", nan_color="#f7f7f7")
    cb4 = fig.colorbar(pcm4, ax=ax, orientation="vertical", shrink=0.8,
                       pad=0.015, aspect=8)
    cb4.set_label("Mean tree distance", fontsize=7)
    cb4.ax.tick_params(labelsize=6)

    # ── Track 5: Gene annotations ─────────────────────────────────────────
    ax = axes[4]
    cat_colors = {
        "MGE_integrase": "#e74c3c",
        "MGE_element": "#f39c12",
        "landmark": "#8e44ad",
    }

    mge = features_df[features_df["category"] != "core"]
    for _, feat in mge.iterrows():
        s_kb = feat["start"] / 1000
        w_kb = max((feat["end"] - feat["start"]) / 1000, genome_kb * 0.002)
        color = cat_colors.get(feat["category"], "#aaa")
        ax.barh(0, w_kb, left=s_kb, height=0.6, color=color,
                edgecolor="none", alpha=0.85)

        gene = feat["gene"]
        if gene and any(kw in gene.lower() for kw in LANDMARK_GENES):
            ax.text(s_kb + w_kb / 2, 0.85, gene, fontsize=5.5,
                    ha="center", va="bottom", rotation=60, style="italic",
                    fontweight="bold")

    ax.set_ylim(-0.3, 1.5)
    ax.set_yticks([])
    ax.set_ylabel("MGE\nfeatures", fontsize=9, rotation=0, ha="right",
                  va="center", labelpad=10)
    ax.set_xlabel("Genome position (kb)", fontsize=10)
    ax.set_xlim(0, genome_kb)
    sns.despine(ax=ax, left=True)

    ax.legend(
        handles=[
            Patch(color="#e74c3c", label="Integrase / transposase"),
            Patch(color="#f39c12", label="MGE element"),
            Patch(color="#8e44ad", label="Landmark gene"),
        ],
        fontsize=6.5, frameon=False, loc="upper right", ncol=3,
    )
    # Assuming ax5 is your bottom track and 'mappable' is the output 
    # from one of your top pcolormesh/imshow plots (e.g., track 1)
    cbar = fig.colorbar(cb4, ax=ax)
    cbar.ax.set_visible(False) # Hide the colorbar, but keep the space it takes up

    # Despine heatmap tracks (they don't need spines)
    for a in axes[:4]:
        sns.despine(ax=a, left=True, bottom=True)
        a.tick_params(axis="x", length=0)

    fig.savefig(output_path, format="svg")
    log.info(f"  Saved {output_path}")
    plt.close()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Map pcDBG boundary severity to reference genome (v3, heatmap)"
    )
    parser.add_argument("--severity", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--reference", required=True, help="GenBank file")
    parser.add_argument("--output", required=True, help="Output prefix")
    parser.add_argument("--window", type=int, default=2000)
    parser.add_argument("--step", type=int, default=500)
    parser.add_argument("--max-ref-dist", type=int, default=10000,
                        help="Max bp between unitig midpoints for local edges")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--evalue", type=float, default=1e-10)
    parser.add_argument("--min-pident", type=float, default=90.0)
    parser.add_argument("--min-coverage", type=float, default=0.5)
    parser.add_argument("--flank", type=int, default=2500)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    # ── Unitigs ───────────────────────────────────────────────────────────
    unitig_fasta = args.output + "_unitigs.fasta"
    if not os.path.exists(unitig_fasta):
        extract_unitig_sequences(args.db, unitig_fasta)
    else:
        log.info(f"Reusing {unitig_fasta}")

    # ── BLAST ─────────────────────────────────────────────────────────────
    ref_fasta = args.output + "_reference.fasta"
    log.info("Parsing reference...")
    record = extract_reference_fasta(args.reference, ref_fasta)
    genome_length = len(record.seq)
    ref_name = record.description[:80] if record.description else record.id

    blast_tsv = args.output + "_blast_hits.tsv"
    if not os.path.exists(blast_tsv):
        run_blast(unitig_fasta, ref_fasta, blast_tsv,
                  threads=args.threads, evalue=args.evalue)
    else:
        log.info(f"Reusing {blast_tsv}")

    log.info("Parsing BLAST hits...")
    hits = parse_blast_hits(blast_tsv, args.min_pident, args.min_coverage)
    log.info(f"  {len(hits):,} unitigs mapped")

    # ── Classify edges ────────────────────────────────────────────────────
    log.info("Classifying edges...")
    lc, ld, struct, unmap = classify_edges(
        args.severity, hits, max_ref_dist=args.max_ref_dist,
    )

    # ── Windows ───────────────────────────────────────────────────────────
    wdf = compute_window_stats(lc, ld, struct, unmap,
                               genome_length, args.window, args.step)
    wdf.to_csv(args.output + "_window_stats.tsv", sep="\t", index=False,
               float_format="%.6f")

    total_local = len(lc) + len(ld)
    log.info(f"  Local: {total_local:,} ({len(lc):,} conc, {len(ld):,} disc)")
    if total_local > 0:
        log.info(f"  Local concordance: {len(lc)/total_local:.4f}")
    log.info(f"  Structural: {len(struct)//2:,}  |  Unmapped-boundary: {len(unmap):,}")

    # ── Features ──────────────────────────────────────────────────────────
    log.info("Parsing features...")
    fdf = parse_features(args.reference)

    # ── Enrichment ────────────────────────────────────────────────────────
    log.info("Computing enrichment...")
    edf = compute_enrichment(wdf, fdf, genome_length, args.step, flank=args.flank)
    if len(edf) > 0:
        edf.to_csv(args.output + "_feature_enrichment.tsv", sep="\t",
                    index=False, float_format="%.4f")
        log.info("  Top 15 by discordance enrichment:")
        for _, r in edf.head(15).iterrows():
            log.info(f"    {r['gene']:15s} {r['product'][:40]:40s} "
                     f"fold={r['fold_enrichment']:.2f} disc={r['discordant_edges']}")

    # ── Plot ──────────────────────────────────────────────────────────────
    plot_genome(wdf, fdf, genome_length, args.output + "_genome_heatmap.svg",
                reference_name=ref_name)

    # ── Unitig map ────────────────────────────────────────────────────────
    with open(args.output + "_unitig_map.tsv", "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["unitig_id", "ref_start", "ref_end", "ref_midpoint"])
        for uid, (rs, re, mid) in sorted(hits.items(), key=lambda x: x[1][2]):
            w.writerow([uid, rs, re, f"{mid:.0f}"])

    log.info("Done.")


if __name__ == "__main__":
    main()
