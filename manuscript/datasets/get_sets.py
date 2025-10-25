#!/usr/bin/env python3
# download_genomes.py
#
# Build 3 x N genome sets with increasing diversity, with caching:
# - Skips re-downloading metadata if already present (unless --force-metadata).
# - Skips re-downloading candidates if target folders already contain enough .fna.gz (unless --force-download).
# - Skips Mash sketch/dist if an up-to-date signature is present (unless --force-sketch, or k/sketch changed).
# Then selects SetL/SetM/SetH and creates seaborn plots.

import argparse, csv, gzip, io, os, re, shutil, subprocess, sys, math, time, hashlib
import numpy as np
import pandas as pd
from collections import defaultdict, Counter
from statistics import median
from random import Random

# ------- Defaults (override via CLI) -------
N = 500
THREADS = 8

# Mash defaults tuned to reduce "all 1.0" saturation at long distances
MASH_K = 17           # was 21
MASH_SKETCH = 100000  # was 10k

WORK = "work"
SETS = "sets"
FIGS = "figs"

PLASMIDS_POLICY = "keep"  # "keep" or "drop" (drop => chromosomes only via header heuristic)

# Distance thresholds / selection targets
L_THRESH_START = 0.0015
L_THRESH_STEP  = 0.0010
L_THRESH_MAX   = 0.0100

M_CLUSTER_CUT  = 0.0100
M_MIN_FLOOR    = 0.0050    # start floor
M_MIN_FLOOR_LO = 0.0020    # how far we may relax the floor
M_MAX_CAP      = 0.0600    # keep distances “within-species”

# Set H spread controls (min floor + max cap)
H_MIN_DIST     = 0.10
H_MAX_DIST     = 0.95

# Candidate caps
LM_MAX_CANDIDATES = 6000
H_MAX_CANDIDATES  = 2000

# Length matching across sets
LENGTH_MATCH = True
LENGTH_TOLERANCE = 0.10  # allow +/-10%

# Optional SetH restriction to a phylum (name or taxid). Use with --taxdump.
H_PHYLUM = None

# Force/caching flags (set via CLI)
FORCE_METADATA = False
FORCE_DOWNLOAD = False
FORCE_SKETCH   = False

RNG = Random(42)
LENS = {}  # filled at runtime

# =========================== Utilities ===========================
def sh(cmd, check=True):
    print(f"[cmd] {cmd}")
    return subprocess.run(cmd, shell=True, check=check)

def ensure_dir(p): os.makedirs(p, exist_ok=True)

def exists_nonempty(p): return os.path.isfile(p) and os.path.getsize(p) > 0

def download(url, out):
    ensure_dir(os.path.dirname(out))
    sh(f"curl -L {url} -o {out}")

def load_assembly_summary(path):
    header = None
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("#") and "assembly_accession" in ln and "organism_name" in ln:
                header = ln.lstrip("#").strip().split("\t")
                break
        if not header:
            raise RuntimeError("Header not found in assembly_summary.")
        for ln in f:
            if ln.startswith("#"):  # skip other comments
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) != len(header):
                continue
            rows.append(dict(zip(header, parts)))
    return rows

def curated_latest_complete(rows):
    out = []
    for r in rows:
        rc = (r.get("refseq_category") or "").lower()
        lvl = r.get("assembly_level") or ""
        vs  = r.get("version_status") or ""
        ftp = r.get("ftp_path") or ""
        if vs != "latest": continue
        if lvl not in {"Complete Genome", "Chromosome"}: continue
        if rc not in {"reference genome", "representative genome"}: continue
        if not (ftp.startswith("ftp://") or ftp.startswith("https://")): continue
        out.append(r)
    return out

def latest_complete_any(rows):
    out = []
    for r in rows:
        vs  = r.get("version_status") or ""
        lvl = r.get("assembly_level") or ""
        ftp = r.get("ftp_path") or ""
        if vs != "latest": continue
        if lvl not in {"Complete Genome", "Chromosome"}: continue
        if not (ftp.startswith("ftp://") or ftp.startswith("https://")): continue
        out.append(r)
    return out

def score_row(r):
    rc = (r.get("refseq_category") or "").lower()
    lvl = r.get("assembly_level") or ""
    typ = 1 if "type" in (r.get("relation_to_type_material") or "").lower() else 0
    date = r.get("seq_rel_date") or "1970-01-01"
    date_key = tuple(int(x) for x in re.split(r"[-/]", date.strip())[:3] if x.isdigit())
    rc_s  = {"reference genome": 3, "representative genome": 2}.get(rc, 1)
    lvl_s = {"Complete Genome": 4, "Chromosome": 3}.get(lvl, 0)
    return (rc_s, lvl_s, typ, date_key)

def genus_name(organism_name):
    parts = (organism_name or "").split()
    if not parts: return ""
    return parts[1] if parts[0].lower() == "candidatus" and len(parts) > 1 else parts[0]

def accession_from_fname(p):
    m = re.search(r"(GC[AF]_\d+\.\d+)", os.path.basename(p))
    return m.group(1) if m else os.path.basename(p)

def fasta_total_len(path):
    opener = gzip.open if path.endswith(".gz") else open
    tot = 0
    with opener(path, "rt", encoding="utf-8", errors="ignore") as f:
        keep = True
        for ln in f:
            if ln.startswith(">"):
                if PLASMIDS_POLICY == "drop" and re.search(r"plasmid", ln, flags=re.I):
                    keep = False
                else:
                    keep = True
                continue
            if keep:
                tot += len(ln.strip())
    return tot

def quantile_bins(lengths, n_bins=10):
    xs = sorted(lengths)
    if not xs: return []
    edges = []
    for q in range(1, n_bins):
        k = int(round(q * len(xs) / n_bins))
        edges.append(xs[max(0, min(k, len(xs)-1))])
    return edges

def bin_index(val, edges):
    b = 0
    for e in edges:
        if val <= e: return b
        b += 1
    return b

def recursive_fna_list(in_dir):
    files = []
    for root, _, fnames in os.walk(in_dir):
        for f in fnames:
            if f.endswith(".fna.gz"):
                files.append(os.path.join(root, f))
    return sorted(files)

# =========================== Taxonomy helpers (optional) ===========================
def parse_taxdump(taxdump_dir):
    """Return maps: parent[taxid], rank[taxid], sci[taxid], name2taxid[name.lower()]."""
    parent = {}
    rank = {}
    with open(os.path.join(taxdump_dir, "nodes.dmp"), "r", encoding="utf-8") as f:
        for ln in f:
            parts = [p.strip() for p in ln.split("|")]
            if len(parts) >= 3:
                tid = int(parts[0]); par = int(parts[1]); rk = parts[2]
                parent[tid] = par; rank[tid] = rk
    sci = {}
    name2taxid = {}
    with open(os.path.join(taxdump_dir, "names.dmp"), "r", encoding="utf-8") as f:
        for ln in f:
            parts = [p.strip() for p in ln.split("|")]
            if len(parts) >= 4:
                tid = int(parts[0]); name = parts[1]; cls = parts[3]
                if cls == "scientific name":
                    sci[tid] = name
                    name2taxid[name.lower()] = tid
    return parent, rank, sci, name2taxid

def taxid_for(name_or_id, name2taxid):
    s = str(name_or_id).strip()
    if re.fullmatch(r"\d+", s):
        return int(s)
    return name2taxid.get(s.lower())

def is_descendant(child_taxid, ancestor_taxid, parent):
    a = int(ancestor_taxid); t = int(child_taxid)
    seen = set()
    while t not in seen and t != 1:
        if t == a: return True
        seen.add(t)
        t = parent.get(t, 1)
    return t == a

# =========================== Mash helpers with caching ===========================
def filelist_signature(filepaths, params_str):
    """
    Build a stable signature over the file list and key params (k, sketch).
    Uses path, size, and mtime to detect changes.
    """
    h = hashlib.sha256()
    h.update(params_str.encode())
    for p in filepaths:
        try:
            st = os.stat(p)
            h.update(p.encode())
            h.update(str(st.st_size).encode())
            h.update(str(int(st.st_mtime)).encode())
        except FileNotFoundError:
            # missing file: still include path to change signature
            h.update(p.encode())
            h.update(b"0"); h.update(b"0")
    return h.hexdigest()

def mash_sketch_and_dist(in_dir, tag, threads, force=False):
    """
    Create/refresh Mash sketch+dist.
    Skips if WORK/{tag}.msh and WORK/{tag}.dist.tsv exist and signature matches current file list + params.
    """
    ensure_dir(WORK)
    msh   = os.path.join(WORK, f"{tag}.msh")
    dist  = os.path.join(WORK, f"{tag}.dist.tsv")
    flist = os.path.join(WORK, f"{tag}.files.txt")
    sigf  = os.path.join(WORK, f"{tag}.sig")

    files = recursive_fna_list(in_dir)
    with open(flist, "w") as out:
        for p in files: out.write(p + "\n")

    params_str = f"k={MASH_K};sketch={MASH_SKETCH}"
    sig_now = filelist_signature(files, params_str)

    if (not force) and exists_nonempty(msh) and exists_nonempty(dist) and exists_nonempty(sigf):
        prev = open(sigf).read().strip()
        if prev == sig_now:
            print(f"[cache] Mash sketch+dist for {tag} up-to-date; skipping.")
            return dist

    # (re)compute
    print(f"[mash] Building sketch for {tag} (k={MASH_K}, sketch={MASH_SKETCH}, n_files={len(files)})")
    sh(f"mash sketch -k {MASH_K} -s {MASH_SKETCH} -p {threads} -o {WORK}/{tag} -l {flist}")
    print(f"[mash] Computing self-dist for {tag}")
    sh(f"mash dist -p {threads} {msh} {msh} > {dist}")
    with open(sigf, "w") as out: out.write(sig_now + "\n")
    return dist

def read_dist(dist_path):
    adj = defaultdict(dict)
    names = set()
    with open(dist_path, "r") as f:
        for ln in f:
            parts = ln.strip().split("\t")
            if len(parts) < 3: 
                continue
            a,b,d = parts[0], parts[1], float(parts[2])
            names.add(a); names.add(b)
            if a != b:
                adj[a][b] = d
                adj[b][a] = d
    return sorted(names), adj

# =========================== Graph helpers ===========================
class DSU:
    def __init__(self, items):
        self.parent = {x:x for x in items}
        self.rank = {x:0 for x in items}
    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x
    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb: return
        if self.rank[ra] < self.rank[rb]:
            self.parent[ra] = rb
        elif self.rank[rb] < self.rank[ra]:
            self.parent[rb] = ra
        else:
            self.parent[rb] = ra
            self.rank[ra] += 1

def components_from_threshold(names, adj, thr):
    dsu = DSU(names)
    for i in names:
        for j, d in adj[i].items():
            if d <= thr:
                dsu.union(i, j)
    comp = defaultdict(list)
    for x in names:
        comp[dsu.find(x)].append(x)
    return list(comp.values())

# =========================== Selectors ===========================
def farthest_first(cands, adj, min_floor, max_cap, quota=None, rng=RNG, target_n=None):
    """
    Maximin selection with both a minimum distance floor and a maximum cap.
    Skips candidates whose current min-distance-to-set is > max_cap to avoid extremes.
    """
    if not cands: return []
    if target_n is None: target_n = N
    seed = rng.choice(cands)
    selected = [seed]

    bin_counts = defaultdict(int)
    getbin = (lambda _: 0)
    if quota:
        getbin = quota["_getbin"]
        bin_counts[getbin(seed)] += 1

    cur_min = {x: (adj[x][seed] if seed in adj[x] else 1.0) for x in cands}
    cur_min[seed] = 0.0

    stalls = 0
    while len(selected) < target_n:
        pool = []
        for x in cands:
            if x in cur_min and x not in selected:
                if quota:
                    b = getbin(x)
                    lim = quota["limits"].get(b, 10**9)
                    if bin_counts[b] >= int((1+LENGTH_TOLERANCE)*lim):
                        continue
                pool.append(x)
        if not pool: break
        nxt = max(pool, key=lambda x: cur_min[x])

        # Cap the extreme candidates
        if cur_min[nxt] > max_cap:
            cur_min.pop(nxt, None)
            stalls += 1
            if stalls > 5 and max_cap < 0.999:
                max_cap = min(0.999, max_cap + 0.02)  # gentle relaxation
            continue

        # Enforce minimum spread
        if cur_min[nxt] < min_floor:
            if quota:
                b = getbin(nxt)
                if bin_counts[b] < int((1+LENGTH_TOLERANCE)*quota["limits"].get(b, 0)):
                    pass
                else:
                    cur_min.pop(nxt, None)
                    continue
            else:
                cur_min.pop(nxt, None)
                continue

        selected.append(nxt)
        if quota:
            bin_counts[getbin(nxt)] += 1

        # Update distances
        for x in cands:
            if x in cur_min and x not in selected:
                dx = adj[x].get(nxt, 1.0)
                if dx < cur_min[x]:
                    cur_min[x] = dx

        if len(selected) % 50 == 0:
            print(f"[farthest] {len(selected)} selected; min(d,new)={cur_min[nxt]:.3f} (cap={max_cap:.2f})")
    return selected[:target_n]

def farthest_first_seeded(cands, adj, seed_set=None, min_floor=0.0, max_cap=1.0,
                          quota=None, target_n=None, rng=RNG):
    """
    Farthest-first (maximin) inside a candidate set, starting from an optional seed_set.
      - Enforces a minimum distance floor (min_floor)
      - Skips extreme outliers whose min-dist-to-set > max_cap
      - Honors optional length-decile quota (softly)
      - Returns up to target_n selections (guaranteed to stop)
    """
    if not cands:
        return []
    if target_n is None:
        target_n = N

    selected = list(seed_set) if seed_set else [rng.choice(cands)]
    suppressed = set()  # candidates we decided to skip (failed floor or cap)

    bin_counts = defaultdict(int)
    getbin = (lambda _: 0)
    if quota:
        getbin = quota["_getbin"]
        for s in selected:
            bin_counts[getbin(s)] += 1

    # current min distance for each candidate to the selected set
    cur_min = {}
    for x in cands:
        if x in selected:
            cur_min[x] = 0.0
        else:
            dmin = 1.0
            for s in selected:
                dmin = min(dmin, adj[x].get(s, 1.0))
            cur_min[x] = dmin

    stalls = 0
    while len(selected) < target_n:
        pool = []
        for x in cands:
            if x in selected or x in suppressed:
                continue
            v = cur_min.get(x, 1.0)
            # hard constraints
            if v > max_cap:
                suppressed.add(x)
                continue
            if quota:
                b = getbin(x)
                lim = quota["limits"].get(b, 10**9)
                if bin_counts[b] >= int((1 + LENGTH_TOLERANCE) * lim):
                    continue
            pool.append(x)

        if not pool:
            break

        nxt = max(pool, key=lambda x: cur_min.get(x, 1.0))
        vnx = cur_min.get(nxt, 1.0)

        # enforce min-floor
        if vnx < min_floor:
            suppressed.add(nxt)     # mark and move on
            stalls += 1
            if stalls > 50:         # give up if spinning
                break
            continue

        # accept nxt
        selected.append(nxt)
        if quota:
            bin_counts[getbin(nxt)] += 1

        # update distances
        for x in cands:
            if x in selected or x in suppressed:
                continue
            dx = adj[x].get(nxt, 1.0)
            if dx < cur_min.get(x, 1.0):
                cur_min[x] = dx

        stalls = 0
        if len(selected) % 50 == 0:
            print(f"[M farthest] {len(selected)} selected; min(d,new)={vnx:.3f} (cap={max_cap:.2f})")

    return selected[:target_n]


def stratified_with_floor(clusters, adj, min_floor, target_n, rng=RNG, quota=None):
    selected = []
    for cl in clusters:
        rng.shuffle(cl)
    idxs = list(range(len(clusters)))

    bin_counts = defaultdict(int)
    getbin = (lambda _: 0)
    if quota:
        getbin = quota["_getbin"]

    while len(selected) < target_n and any(cl for cl in clusters):
        progressed = False
        for ci in idxs:
            cl = clusters[ci]
            while cl:
                g = cl.pop()
                if quota:
                    b = getbin(g)
                    lim = quota["limits"].get(b, 10**9)
                    if bin_counts[b] >= int((1+LENGTH_TOLERANCE)*lim):
                        continue
                ok = True
                for h in selected:
                    d = adj[g].get(h, 1.0)
                    if d < min_floor:
                        ok = False; break
                if ok:
                    selected.append(g)
                    if quota:
                        bin_counts[getbin(g)] += 1
                    progressed = True
                    break
            if len(selected) >= target_n:
                break
        if not progressed:
            break
    return selected[:target_n]

# =========================== Pipeline steps (with caching) ===========================
def step_fetch_metadata():
    ensure_dir(WORK)
    out = os.path.join(WORK, "assembly_summary_refseq_bacteria.txt")
    if exists_nonempty(out) and not FORCE_METADATA:
        print("[cache] assembly_summary present; skipping download.")
        return out
    download("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt", out)
    return out

def step_pick_candidates(asm_path, taxmaps=None, h_phylum=None):
    """
    L/M: choose top species by count from all high-quality (latest+complete/chromosome) RefSeq.
    H: curated one-per-genus; optionally restrict to a phylum (name or taxid) if taxmaps provided.
    """
    rows  = load_assembly_summary(asm_path)
    cur   = curated_latest_complete(rows)   # for H (curated)
    hqany = latest_complete_any(rows)       # for L/M (broad within-species)

    # L/M: species with most high-quality assemblies
    by_species = defaultdict(list)
    for r in hqany:
        sp = r.get("species_taxid") or ""
        if sp: by_species[sp].append(r)

    species_ranked = sorted(by_species.items(), key=lambda kv: len(kv[1]), reverse=True)
    need = max(N*3, 1000)
    top_sp, top_list = next(((sp, lst) for sp,lst in species_ranked if len(lst) >= need), species_ranked[0])

    sp_name = " ".join((top_list[0]["organism_name"] or "").split()[:2])
    print(f"[candidates] Chosen species for L/M: {sp_name} (species_taxid={top_sp}, {len(top_list)} assemblies).")

    RNG.shuffle(top_list)
    if len(top_list) > LM_MAX_CANDIDATES:
        top_list = top_list[:LM_MAX_CANDIDATES]
    lm_accessions = [r["assembly_accession"] for r in top_list]

    # H: curated one-per-genus; optional phylum restriction
    if h_phylum and taxmaps:
        parent, rank, sci, name2taxid = taxmaps
        phyl_id = taxid_for(h_phylum, name2taxid)
        if phyl_id:
            print(f"[candidates] Restricting H pool to phylum '{h_phylum}' (taxid={phyl_id})")
            cur = [r for r in cur
                   if r.get("species_taxid") and is_descendant(int(r["species_taxid"]), phyl_id, parent)]
        else:
            print(f"[warn] Could not resolve phylum '{h_phylum}'. Proceeding without restriction.")

    best_per_species = {}
    for r in cur:
        sp = r.get("species_taxid")
        if not sp: continue
        if sp not in best_per_species or score_row(r) > score_row(best_per_species[sp]):
            best_per_species[sp] = r
    best_per_genus = {}
    for r in best_per_species.values():
        g = genus_name(r["organism_name"])
        if not g: continue
        if g not in best_per_genus or score_row(r) > score_row(best_per_genus[g]):
            best_per_genus[g] = r
    genus_list = list(best_per_genus.values())
    print(f"[candidates] One-per-genus curated H-pool size: {len(genus_list)}")
    RNG.shuffle(genus_list)
    if len(genus_list) > H_MAX_CANDIDATES:
        genus_list = genus_list[:H_MAX_CANDIDATES]
    h_accessions = [r["assembly_accession"] for r in genus_list]

    ensure_dir(os.path.join(WORK, "candidates"))
    lm_list = os.path.join(WORK, "candidates", "LM_accessions.txt")
    h_list  = os.path.join(WORK, "candidates", "H_accessions.txt")
    with open(lm_list, "w") as f: f.write("\n".join(lm_accessions) + "\n")
    with open(h_list,  "w") as f: f.write("\n".join(h_accessions) + "\n")

    with open(os.path.join(WORK, "candidates", "LM_species.txt"), "w") as f:
        f.write(f"{sp_name}\n")
    return lm_list, h_list

def count_fnas(in_dir):
    return len(recursive_fna_list(in_dir))

def num_lines(path):
    return sum(1 for _ in open(path, "r") if _.strip())

def step_download_candidates(lm_list, h_list):
    out_lm = os.path.join(WORK, "candidates_LM")
    out_h  = os.path.join(WORK, "candidates_H")
    ensure_dir(out_lm); ensure_dir(out_h)

    # Skip logic: if enough files already present and not forcing, skip the NGD call.
    lm_need = num_lines(lm_list)
    h_need  = num_lines(h_list)

    have_lm = count_fnas(out_lm)
    have_h  = count_fnas(out_h)

    base_lm = "--section refseq --formats fasta --assembly-levels complete,chromosome --flat-output"
    base_h  = "--section refseq --formats fasta --assembly-levels complete,chromosome --refseq-categories reference,representative --flat-output"

    if (have_lm >= lm_need) and not FORCE_DOWNLOAD:
        print(f"[cache] L/M candidates present ({have_lm}/{lm_need}); skipping download.")
    else:
        sh(f"ncbi-genome-download {base_lm} --parallel {THREADS} -A {lm_list} --output-folder {out_lm} bacteria")

    if (have_h >= h_need) and not FORCE_DOWNLOAD:
        print(f"[cache] H candidates present ({have_h}/{h_need}); skipping download.")
    else:
        sh(f"ncbi-genome-download {base_h}  --parallel {THREADS} -A {h_list}  --output-folder {out_h}  bacteria")

    return out_lm, out_h

def step_lengths(in_dir):
    files = recursive_fna_list(in_dir)
    return {p: fasta_total_len(p) for p in files}

def quotas_from_deciles(lengths_selected_L):
    edges = quantile_bins(list(lengths_selected_L.values()), 10)
    lim = Counter([bin_index(v, edges) for v in lengths_selected_L.values()])
    return edges, lim

def build_quota_struct(edges, lim, target_n=N):
    def _getbin(name_or_path):
        return bin_index(LENS[name_or_path], edges)
    return {"limits": lim, "_getbin": _getbin, "_target_n": target_n}

def stage_set(names, tag):
    ensure_dir(os.path.join(SETS, tag))
    accs = [accession_from_fname(x) for x in names]
    with open(os.path.join(SETS, tag, "accessions.txt"), "w") as f:
        f.write("\n".join(accs) + "\n")
    for p in names:
        dst = os.path.join(SETS, tag, os.path.basename(p))
        if os.path.exists(dst): 
            continue
        try:
            os.symlink(os.path.abspath(p), dst)
        except Exception:
            shutil.copy2(p, dst)

def step_select_sets(lm_dir, h_dir):
    global LENS
    ensure_dir(SETS)

    # lengths (for quotas)
    LENS = {}
    LENS.update(step_lengths(lm_dir))
    LENS.update(step_lengths(h_dir))

    # LM distances (cached)
    lm_dist = mash_sketch_and_dist(lm_dir, "LM", THREADS, force=FORCE_SKETCH)
    names_lm, adj_lm = read_dist(lm_dist)

    # ---- Set L: auto-tune dense threshold ----
    L_names = []
    thr = L_THRESH_START
    while thr <= L_THRESH_MAX:
        comps = components_from_threshold(names_lm, adj_lm, thr)
        big = max(comps, key=len) if comps else []
        if len(big) >= N:
            L_names = RNG.sample(big, N)
            print(f"[L] threshold={thr:.4f} produced component size={len(big)}. Selected {N}.")
            break
        thr += L_THRESH_STEP
    if not L_names:
        raise RuntimeError("Could not find a dense enough cluster for Set L; increase L_THRESH_MAX or LM_MAX_CANDIDATES.")

    # build quota from L lengths (optional)
    edges_L, lim_L = quotas_from_deciles({p:LENS[p] for p in L_names})
    quota_M = build_quota_struct(edges_L, lim_L, N) if LENGTH_MATCH else None
    quota_H = build_quota_struct(edges_L, lim_L, N) if LENGTH_MATCH else None

    # ---- Set M (within-species, diverse) ----
    # Stage A: farthest-first with target floor + within-species cap
    M_names = farthest_first_seeded(
        names_lm, adj_lm,
        seed_set=None,
        min_floor=M_MIN_FLOOR,
        max_cap=M_MAX_CAP,
        quota=quota_M,
        target_n=N,
        rng=RNG
    )

    # Stage B: if we stalled, gradually relax the floor (still keep the cap)
    if len(M_names) < N:
        print(f"[M] Only {len(M_names)} with floor={M_MIN_FLOOR:.4f}; relaxing floor...")
        floor = max(M_MIN_FLOOR * 0.8, M_MIN_FLOOR_LO)
        while len(M_names) < N and floor >= M_MIN_FLOOR_LO - 1e-9:
            M_names = farthest_first_seeded(
                [x for x in names_lm if x not in M_names], adj_lm,
                seed_set=M_names,
                min_floor=floor,
                max_cap=M_MAX_CAP,
                quota=quota_M,
                target_n=N,
                rng=RNG
            )
            floor *= 0.9

    # Stage C: guaranteed fill (drop quota & floor, but still within species cap)
    if len(M_names) < N:
        print(f"[M] Filling the last {N - len(M_names)} without floor/quota (cap stays {M_MAX_CAP:.3f}).")
        pool = [x for x in names_lm if x not in M_names]
        RNG.shuffle(pool)
        for x in pool:
            if len(M_names) >= N:
                break
            # keep cap so we don't slip into cross-species-like extremes
            dmin = min(adj_lm[x].get(y, 1.0) for y in M_names) if M_names else 1.0
            if dmin <= M_MAX_CAP:
                M_names.append(x)


    # ---- H distances & selection (cached + spread caps) ----
    h_dist = mash_sketch_and_dist(h_dir, "H", THREADS, force=FORCE_SKETCH)
    names_h, adj_h = read_dist(h_dist)
    H_names = farthest_first(names_h, adj_h, H_MIN_DIST, H_MAX_DIST, quota=quota_H, rng=RNG, target_n=N)
    if len(H_names) < N:
        print(f"[H] Only got {len(H_names)} with constraints (min={H_MIN_DIST}, max={H_MAX_DIST}); filling with relaxed floor.")
        pool = [x for x in names_h if x not in H_names]
        RNG.shuffle(pool)
        for x in pool:
            if len(H_names) >= N: break
            ok = all(adj_h[x].get(y,1.0) >= max(0.5*H_MIN_DIST, 0.05) for y in H_names)
            if ok: H_names.append(x)
        H_names = H_names[:N]

    # ---- stage outputs ----
    stage_set(L_names, "SetL")
    stage_set(M_names, "SetM")
    stage_set(H_names, "SetH")

    # ---- summaries ----
    def pairwise_stats(names, adj):
        ds = []
        for i,a in enumerate(names):
            for b in names[i+1:]:
                if b in adj[a]:
                    ds.append(adj[a][b])
        med = median(ds) if ds else float("nan")
        q1 = sorted(ds)[len(ds)//4] if ds else float("nan")
        q3 = sorted(ds)[3*len(ds)//4] if ds else float("nan")
        total_bp = sum(LENS[x] for x in names)
        med_len = median([LENS[x] for x in names])
        return med, q1, q3, total_bp, med_len, len(ds)

    rows = []
    for tag, names, adj in [("SetL", L_names, adj_lm), ("SetM", M_names, adj_lm), ("SetH", H_names, adj_h)]:
        med, q1, q3, total_bp, med_len, npairs = pairwise_stats(names, adj)
        rows.append([tag, len(names), total_bp, med_len, med, q1, q3, npairs])

    with open(os.path.join(SETS, "summary.tsv"), "w") as f:
        f.write("set\tn_genomes\ttotal_bp\tmedian_len\tmedian_mash\tq1_mash\tq3_mash\tn_pairs\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")

    print("[done] Wrote sets/Set*/accessions.txt, staged FASTAs, and sets/summary.tsv")

# =========================== Plotting (seaborn) ===========================
def list_set_basenames(set_dir):
    return {f for f in os.listdir(set_dir) if f.endswith(".fna.gz")}

def read_mash_pairs_for_set(dist_path, set_basenames, max_pairs=None):
    dists = []
    seen = set()
    with open(dist_path, "r") as f:
        for ln in f:
            parts = ln.strip().split("\t")
            if len(parts) < 3: 
                continue
            a, b = parts[0], parts[1]
            try:
                d = float(parts[2])
            except:
                continue
            ba, bb = os.path.basename(a), os.path.basename(b)
            if ba == bb: 
                continue
            if (ba in set_basenames) and (bb in set_basenames):
                key = (ba, bb) if ba < bb else (bb, ba)
                if key in seen: 
                    continue
                seen.add(key)
                dists.append(d)
                if max_pairs and len(dists) >= max_pairs:
                    break
    return dists

def parse_assembly_summary_df(path):
    header = None
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("#") and "assembly_accession" in ln and "organism_name" in ln:
                header = ln.lstrip("#").strip().split("\t")
                break
        if header is None:
            raise RuntimeError("Header not found in assembly_summary.")
        for ln in f:
            if ln.startswith("#"): 
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) != len(header):
                continue
            rows.append(dict(zip(header, parts)))
    import pandas as pd
    df = pd.DataFrame(rows)
    def genus_from_name(name):
        parts = (name or "").split()
        if not parts: return ""
        return parts[1] if parts[0].lower()=="candidatus" and len(parts)>1 else parts[0]
    df["genus"] = df["organism_name"].apply(genus_from_name)
    return df[["assembly_accession","organism_name","species_taxid","genus"]]

def _asm_df_for_tax(work_dir):
    """Parse assembly_summary to the columns we need."""
    path = os.path.join(work_dir, "assembly_summary_refseq_bacteria.txt")
    header = None
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("#") and "assembly_accession" in ln and "organism_name" in ln:
                header = ln.lstrip("#").strip().split("\t")
                break
        if header is None:
            raise RuntimeError("Header not found in assembly_summary.")
        for ln in f:
            if ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) != len(header):
                continue
            rows.append(dict(zip(header, parts)))
    df = pd.DataFrame(rows)
    # Normalize types
    if "species_taxid" in df.columns:
        df["species_taxid"] = pd.to_numeric(df["species_taxid"], errors="coerce").astype("Int64")
    # Genus (from organism_name)
    def _genus(name):
        parts = (name or "").split()
        if not parts: return ""
        return parts[1] if parts[0].lower()=="candidatus" and len(parts)>1 else parts[0]
    df["genus_name"] = df["organism_name"].apply(_genus)
    # Return the useful subset (keep infraspecific_name in case you want strain later)
    keep = ["assembly_accession","organism_name","infraspecific_name","species_taxid","genus_name"]
    return df[keep]

def _ancestor_taxid(taxid, want_rank, parent, rank):
    """Walk up to the wanted rank; return taxid or None."""
    if pd.isna(taxid): return None
    try: t = int(taxid)
    except Exception: return None
    seen = set()
    while t not in seen and t != 1:
        if rank.get(t) == want_rank:
            return t
        seen.add(t)
        t = parent.get(t, 1)
    return None

def _rank_names(species_taxid, parent, rank, sci):
    """Return a dict of names+ids for genus..domain given a species_taxid."""
    out = {
        "species_taxid": species_taxid,
        "species_name": sci.get(int(species_taxid), None) if pd.notna(species_taxid) else None,
        "genus_taxid": None, "genus_name": None,
        "family_taxid": None, "family_name": None,
        "order_taxid": None, "order_name": None,
        "class_taxid": None, "class_name": None,
        "phylum_taxid": None, "phylum_name": None,
        "domain_taxid": None, "domain_name": None,  # NCBI rank: superkingdom
    }
    if pd.isna(species_taxid):
        return out
    # genus → phylum → superkingdom (domain)
    g = _ancestor_taxid(species_taxid, "genus", parent, rank)
    f = _ancestor_taxid(species_taxid, "family", parent, rank)
    o = _ancestor_taxid(species_taxid, "order", parent, rank)
    c = _ancestor_taxid(species_taxid, "class", parent, rank)
    p = _ancestor_taxid(species_taxid, "phylum", parent, rank)
    d = _ancestor_taxid(species_taxid, "superkingdom", parent, rank)

    out.update({
        "genus_taxid": g, "genus_name": sci.get(g),
        "family_taxid": f, "family_name": sci.get(f),
        "order_taxid": o, "order_name": sci.get(o),
        "class_taxid": c, "class_name": sci.get(c),
        "phylum_taxid": p, "phylum_name": sci.get(p),
        "domain_taxid": d, "domain_name": sci.get(d),
    })
    return out

def export_df_tax(sets_dir="sets", work_dir="work", taxdump_dir=None,
                  out_tsv=None, out_parquet=None):
    """
    Build df_tax across SetL/SetM/SetH with rich taxonomy:
      set, assembly_accession, organism_name, infraspecific_name,
      species_taxid, species_name, genus_name (+taxid),
      family/order/class/phylum/domain (names + taxids).
    Writes sets/df_tax.tsv (and .parquet).
    """
    if out_tsv is None:
        out_tsv = os.path.join(sets_dir, "df_tax.tsv")
    if out_parquet is None:
        out_parquet = os.path.join(sets_dir, "df_tax.parquet")

    # 1) Assembly metadata
    asm = _asm_df_for_tax(work_dir)

    # 2) Read accessions per set
    def _acc(tag):
        p = os.path.join(sets_dir, tag, "accessions.txt")
        return [ln.strip() for ln in open(p) if ln.strip()]
    rows = []
    for tag in ["SetL","SetM","SetH"]:
        accs = _acc(tag)
        sub = asm[asm["assembly_accession"].isin(accs)].copy()
        sub.insert(0, "set", tag)
        rows.append(sub)
    df = pd.concat(rows, ignore_index=True)

    # 3) If taxdump present, enrich ranks; else fill species from organism_name fallback
    if taxdump_dir and os.path.isdir(taxdump_dir):
        parent, rank, sci, name2taxid = parse_taxdump(taxdump_dir)
        # Vectorize via .apply (fast enough for 1500 rows)
        taxrecs = df["species_taxid"].apply(lambda tid: _rank_names(tid, parent, rank, sci))
        tax = pd.DataFrame(list(taxrecs.values))
        df_tax = pd.concat([df.reset_index(drop=True), tax.reset_index(drop=True)], axis=1)
        # For rows with missing species_name, fallback to first two tokens of organism_name
        mask = df_tax["species_name"].isna()
        df_tax.loc[mask, "species_name"] = df_tax.loc[mask, "organism_name"].str.split().str[:2].str.join(" ")
    else:
        # No taxdump: provide genus and a naive species string
        df_tax = df.copy()
        df_tax["species_name"] = df_tax["organism_name"].str.split().str[:2].str.join(" ")
        # Fill placeholders
        for col in ["genus_taxid","family_taxid","order_taxid","class_taxid","phylum_taxid","domain_taxid",
                    "family_name","order_name","class_name","phylum_name","domain_name"]:
            df_tax[col] = pd.Series([np.nan]*len(df_tax))

    # 4) Order / save
    cols = [
        "set",
        "assembly_accession",
        "organism_name",
        "infraspecific_name",
        "species_taxid","species_name",
        "genus_taxid","genus_name",
        "family_taxid","family_name",
        "order_taxid","order_name",
        "class_taxid","class_name",
        "phylum_taxid","phylum_name",
        "domain_taxid","domain_name",
    ]
    df_tax = df_tax[cols]
    df_tax.to_csv(out_tsv, sep="\t", index=False)
    try:
        df_tax.to_parquet(out_parquet, index=False)
    except Exception:
        pass  # parquet optional

    print(f"[taxonomy] Wrote {out_tsv}")
    if os.path.exists(out_parquet):
        print(f"[taxonomy] Wrote {out_parquet}")

def shannon(counts):
    tot = sum(counts.values())
    if tot == 0: return 0.0
    H = 0.0
    for c in counts.values():
        p = c / tot
        if p > 0:
            H -= p * math.log(p)
    return H

def make_plots(taxdump=None, max_pairs=None):
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt

    os.makedirs(FIGS, exist_ok=True)
    df_sum = pd.read_csv(os.path.join(SETS, "summary.tsv"), sep="\t")
    print("\n[summary]\n", df_sum, "\n")

    # Distances
    L_b = list_set_basenames(os.path.join(SETS, "SetL"))
    M_b = list_set_basenames(os.path.join(SETS, "SetM"))
    H_b = list_set_basenames(os.path.join(SETS, "SetH"))
    lm_dist_path = os.path.join(WORK, "LM.dist.tsv")
    h_dist_path  = os.path.join(WORK, "H.dist.tsv")

    L_d = read_mash_pairs_for_set(lm_dist_path, L_b, max_pairs=max_pairs)
    M_d = read_mash_pairs_for_set(lm_dist_path, M_b, max_pairs=max_pairs)
    H_d = read_mash_pairs_for_set(h_dist_path,  H_b, max_pairs=max_pairs)

    df_dist = pd.concat([
        pd.DataFrame({"set":"SetL","distance":L_d}),
        pd.DataFrame({"set":"SetM","distance":M_d}),
        pd.DataFrame({"set":"SetH","distance":H_d}),
    ], ignore_index=True)
    print(df_dist.groupby("set")["distance"].describe(), "\n")

    sns.set(context="talk", style="whitegrid")

    g = sns.displot(
        df_dist, x="distance", hue="set", kind="hist",
        element="step", stat="density", common_norm=False, bins=60, height=5, aspect=1.6
    )
    g.set(xlim=(0,1.0), xlabel="Mash distance", ylabel="Density")
    g.fig.tight_layout()
    g.fig.savefig(os.path.join(FIGS, "mash_hist_by_set.png"), dpi=200)

    plt.figure(figsize=(8,5))
    ax = sns.violinplot(data=df_dist, x="set", y="distance", cut=0, inner="quartile")
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("")
    ax.set_ylabel("Mash distance")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "mash_violin_by_set.png"), dpi=200)

    # Taxonomy
    df_asm = parse_assembly_summary_df(os.path.join(WORK, "assembly_summary_refseq_bacteria.txt"))
    sets_df = []
    for tag in ["SetL","SetM","SetH"]:
        accs = [ln.strip() for ln in open(os.path.join(SETS, tag, "accessions.txt")) if ln.strip()]
        sub = df_asm[df_asm["assembly_accession"].isin(accs)].copy()
        sub["set"] = tag
        sets_df.append(sub[["set","assembly_accession","organism_name","species_taxid","genus"]])
    df_tax = pd.concat(sets_df, ignore_index=True)

    # Optional phylum mapping
    if taxdump and os.path.isdir(taxdump):
        try:
            parent, rank, sci, name2taxid = parse_taxdump(taxdump)
            def to_phylum_name(tid):
                try:
                    t = int(str(tid))
                except:
                    return None
                seen = set()
                while t not in seen and t != 1:
                    if rank.get(t) == "phylum":
                        return sci.get(t, None)
                    seen.add(t)
                    t = parent.get(t, 1)
                return None
            df_tax["phylum_name"] = df_tax["species_taxid"].apply(to_phylum_name)
        except Exception as e:
            print("[warn] taxdump parsing failed:", e)

    # Top-N genera per set
    import seaborn as sns
    top_n = 20
    top = (
        df_tax.groupby(["set","genus"])["assembly_accession"]
        .count()
        .reset_index(name="n")
        .sort_values(["set","n"], ascending=[True, False])
        .groupby("set")
        .head(top_n)
    )
    plt.figure(figsize=(12,8))
    ax = sns.barplot(data=top, y="genus", x="n", hue="set", dodge=True)
    ax.set_ylabel(f"Genus (top {top_n})")
    ax.set_xlabel("Count")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, f"taxonomy_top{top_n}_genera.png"), dpi=200)

    # Shannon (genus) – vectorized, no .apply on grouped df
    counts = df_tax.groupby(["set","genus"]).size().unstack(fill_value=0)
    props  = counts.div(counts.sum(axis=1), axis=0).replace(0, np.nan)
    shan_series = -(props * np.log(props)).sum(axis=1).fillna(0.0)
    shan = shan_series.reset_index()
    shan.columns = ["set", "shannon_genus"]

    plt.figure(figsize=(6,4))
    ax = sns.barplot(data=shan, x="set", y="shannon_genus")
    ax.set_ylabel("Shannon diversity (genus)")
    ax.set_xlabel("")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "shannon_genus.png"), dpi=200)

    # Phylum stacked (if available)
    if "phylum_name" in df_tax.columns:
        counts = (
            df_tax.groupby(["set","phylum_name"])["assembly_accession"]
            .count().rename("n").reset_index()
        )
        totals = counts.groupby("set")["n"].sum().rename("tot")
        counts = counts.merge(totals, on="set")
        counts["prop"] = counts["n"] / counts["tot"]
        top_k = 12
        top_phyla = (
            counts.groupby("phylum_name")["n"].sum()
            .sort_values(ascending=False).head(top_k).index.tolist()
        )
        counts["phylum_plot"] = counts["phylum_name"].where(counts["phylum_name"].isin(top_phyla), "Other")
        piv = counts.groupby(["set","phylum_plot"])["prop"].sum().reset_index()
        piv = piv.pivot(index="set", columns="phylum_plot", values="prop").fillna(0.0)
        piv = piv[sorted(piv.columns)]

        ax = piv.plot(kind="bar", stacked=True, figsize=(10,5))
        ax.set_ylabel("Proportion")
        ax.set_xlabel("")
        ax.legend(title="Phylum", bbox_to_anchor=(1.02,1), loc="upper left", borderaxespad=0.)
        plt.tight_layout()
        plt.savefig(os.path.join(FIGS, "phylum_stacked.png"), dpi=200)

    print("\n[wrote]")
    print("  figs/mash_hist_by_set.png")
    print("  figs/mash_violin_by_set.png")
    print("  figs/taxonomy_top20_genera.png")
    print("  figs/shannon_genus.png")
    if "phylum_name" in df_tax.columns:
        print("  figs/phylum_stacked.png")

# =========================== CLI ===========================
def main():
    global N, THREADS, PLASMIDS_POLICY, LENGTH_MATCH, MASH_K, MASH_SKETCH
    global H_MIN_DIST, H_MAX_DIST, H_PHYLUM, FORCE_METADATA, FORCE_DOWNLOAD, FORCE_SKETCH

    ap = argparse.ArgumentParser(description="Build 3xN genome sets (increasing diversity) with caching, then plot.")
    ap.add_argument("--n", type=int, default=N)
    ap.add_argument("--threads", type=int, default=THREADS)
    ap.add_argument("--plasmids", choices=["keep","drop"], default=PLASMIDS_POLICY)
    ap.add_argument("--length-match", action="store_true", default=LENGTH_MATCH,
                    help="Match Set M/H genome-length deciles to Set L")
    ap.add_argument("--mash-k", type=int, default=MASH_K, help="Mash k-mer size (default 17)")
    ap.add_argument("--mash-sketch", type=int, default=MASH_SKETCH, help="Mash sketch size (default 100000)")
    # Set H controls
    ap.add_argument("--h-min-dist", type=float, default=H_MIN_DIST, help="SetH minimum distance floor (default 0.10)")
    ap.add_argument("--h-max-dist", type=float, default=H_MAX_DIST, help="SetH maximum min-distance cap (default 0.95)")
    ap.add_argument("--h-phylum", default=None, help="Restrict SetH candidates to this phylum (name or taxid). Requires --taxdump.")
    # Taxdump for optional phylum restriction & phylum plots
    ap.add_argument("--taxdump", default=None, help="Path to NCBI taxdump dir (nodes.dmp, names.dmp)")
    # Plotting
    ap.add_argument("--no-plots", action="store_true", help="Skip seaborn plots")
    ap.add_argument("--max-pairs", type=int, default=None, help="Cap #pairs read for plots (speed up)")
    # Caching/forcing
    ap.add_argument("--force-metadata", action="store_true", help="Redownload assembly_summary")
    ap.add_argument("--force-download", action="store_true", help="Redownload candidate FASTAs")
    ap.add_argument("--force-sketch", action="store_true", help="Recompute Mash sketches/distances")

    args = ap.parse_args()

    # Set globals from CLI
    N = args.n
    THREADS = args.threads
    PLASMIDS_POLICY = args.plasmids
    LENGTH_MATCH = args.length_match
    MASH_K = args.mash_k
    MASH_SKETCH = args.mash_sketch
    H_MIN_DIST = args.h_min_dist
    H_MAX_DIST = args.h_max_dist
    H_PHYLUM = args.h_phylum
    FORCE_METADATA = args.force_metadata
    FORCE_DOWNLOAD = args.force_download
    FORCE_SKETCH = args.force_sketch

    taxmaps = None
    if args.taxdump:
        try:
            taxmaps = parse_taxdump(args.taxdump)
        except Exception as e:
            print("[warn] Could not parse taxdump:", e)

    ensure_dir(WORK); ensure_dir(SETS); ensure_dir(FIGS)

    asm = step_fetch_metadata()
    lm_list, h_list = step_pick_candidates(asm, taxmaps=taxmaps, h_phylum=H_PHYLUM)
    lm_dir, h_dir = step_download_candidates(lm_list, h_list)
    step_select_sets(lm_dir, h_dir)
    export_df_tax(sets_dir=SETS, work_dir=WORK, taxdump_dir=args.taxdump)


if __name__ == "__main__":
    main()

