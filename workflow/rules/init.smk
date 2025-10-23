from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, Iterable, Optional


configfile: "config.yaml"


#########################################
# Configuration helpers and directories #
#########################################


INPUT_DIR = Path(config["input_dir"]).resolve()
INTERMEDIATE_DIR = Path(config["intermediate_dir"]).resolve()
OUTPUT_DIR = Path(config["output_dir"]).resolve()
KMER_LENGTH = int(config["kmer_length"])
ALLOWED_FASTA_SUFFIXES = {"fa", "fasta", "fna", "ffa"}


def dir_input() -> Path:
    return INPUT_DIR


def dir_intermediate() -> Path:
    return INTERMEDIATE_DIR


def dir_output() -> Path:
    return OUTPUT_DIR


#############################################
# Batch manifest parsing and cached lookups #
#############################################


def _sample_name_from_path(entry: str) -> str:
    """Infer the logical sample name from a FASTA path."""
    filename = Path(entry).name
    if filename.endswith(".gz"):
        filename = filename[:-3]
    sample, dot, suffix = filename.rpartition(".")
    if not dot:
        raise ValueError(f"Unable to determine suffix for '{entry}'.")
    if suffix not in ALLOWED_FASTA_SUFFIXES:
        raise ValueError(
            f"Unknown FASTA suffix '{suffix}' in '{entry}'. "
            f"Expected one of {sorted(ALLOWED_FASTA_SUFFIXES)}."
        )
    return sample


def _load_batches() -> Dict[str, Dict[str, str]]:
    """Parse manifests in the input directory into a batch → sample map."""

    manifests = sorted(INPUT_DIR.glob("*.txt"))
    if not manifests:
        raise FileNotFoundError(
            f"No manifests found in '{INPUT_DIR}'. Provide at least one '*.txt' file."
        )

    batches: Dict[str, Dict[str, str]] = {}
    for manifest in manifests:
        batch = manifest.stem
        tree_path = INPUT_DIR / f"{batch}.nwk"
        if not tree_path.exists():
            raise FileNotFoundError(
                f"Missing Newick tree for batch '{batch}'. Expected '{tree_path}'."
            )

        samples: Dict[str, str] = {}
        with manifest.open() as handle:
            for line_no, raw in enumerate(handle, start=1):
                entry = raw.strip()
                if not entry:
                    continue
                sample_name = _sample_name_from_path(entry)
                samples[sample_name] = entry

        if not samples:
            raise ValueError(f"Manifest '{manifest}' contains no sample entries.")

        batches[batch] = samples

    return batches


BATCH_SAMPLE_MAP = _load_batches()


def get_batches() -> Iterable[str]:
    return tuple(sorted(BATCH_SAMPLE_MAP))


#####################################
# Global files for individual batches
#####################################


def _cuttlefish_root(batch: str) -> Path:
    return INTERMEDIATE_DIR / "cuttlefish" / batch


def _minimal_cuts_root(batch: str) -> Path:
    return INTERMEDIATE_DIR / "minimalcuts"


def _tree_root() -> Path:
    return INTERMEDIATE_DIR / "tree"


def _bubble_root() -> Path:
    return INTERMEDIATE_DIR / "bubblegun"


def _require_batch(batch: Optional[str] = None, _batch: Optional[str] = None) -> str:
    value = batch if batch is not None else _batch
    if value is None:
        raise ValueError("Batch identifier is required.")
    return value


def _require_ext(ext: Optional[str] = None, _ext: Optional[str] = None) -> str:
    value = ext if ext is not None else _ext
    if value is None:
        raise ValueError("File extension is required.")
    return value


def fn_cuttlefish_out(*, batch: Optional[str] = None, _batch: Optional[str] = None, ext: Optional[str] = None, _ext: Optional[str] = None):
    batch_id = _require_batch(batch, _batch)
    extension = _require_ext(ext, _ext)
    return str(_cuttlefish_root(batch_id) / f"{batch_id}_compcoloreddbg_k{KMER_LENGTH}.{extension}")


def fn_colormtx(batch: Optional[str] = None, _batch: Optional[str] = None) -> str:
    batch_id = _require_batch(batch, _batch)
    return str(INTERMEDIATE_DIR / "cuttlefish" / f"{batch_id}_unitigcolors_k{KMER_LENGTH}.tsv")


def fn_uqcolors(batch: Optional[str] = None, _batch: Optional[str] = None) -> str:
    batch_id = _require_batch(batch, _batch)
    return str(INTERMEDIATE_DIR / "cuttlefish" / f"{batch_id}_uqcolors_k{KMER_LENGTH}.tsv")


def fn_redundantcolors(batch: Optional[str] = None, _batch: Optional[str] = None) -> str:
    batch_id = _require_batch(batch, _batch)
    return str(INTERMEDIATE_DIR / "cuttlefish" / f"{batch_id}_redundantcolors_k{KMER_LENGTH}.csv")


def fn_minimalcuts(batch: Optional[str] = None, _batch: Optional[str] = None) -> str:
    batch_id = _require_batch(batch, _batch)
    return str(_minimal_cuts_root(batch_id) / f"{batch_id}_minimalcuts_k{KMER_LENGTH}")


def fn_minimalcuts_plotdir(batch: Optional[str] = None, _batch: Optional[str] = None) -> str:
    batch_id = _require_batch(batch, _batch)
    return str(_minimal_cuts_root(batch_id) / f"{batch_id}_plots_k{KMER_LENGTH}")


def fn_unitig2cuts(batch: Optional[str] = None, _batch: Optional[str] = None) -> str:
    batch_id = _require_batch(batch, _batch)
    return str(_minimal_cuts_root(batch_id) / f"{batch_id}_unitig2cuts_k{KMER_LENGTH}")


def fn_sqldb(batch: Optional[str] = None, _batch: Optional[str] = None) -> str:
    batch_id = _require_batch(batch, _batch)
    return str(INTERMEDIATE_DIR / "sqldb" / f"{batch_id}_k{KMER_LENGTH}.sqldb")


def fn_tree_clean(batch: Optional[str] = None, _batch: Optional[str] = None) -> str:
    batch_id = _require_batch(batch, _batch)
    return str(_tree_root() / f"{batch_id}.nwk")


def fn_tree_sorted(batch: Optional[str] = None, _batch: Optional[str] = None) -> str:
    return fn_tree_clean(batch=batch, _batch=_batch)


def fn_tree_dirty(batch: Optional[str] = None, _batch: Optional[str] = None) -> str:
    batch_id = _require_batch(batch, _batch)
    return str(_tree_root() / f"{batch_id}.nwk_dirty")


def fn_leaves_sorted(batch: Optional[str] = None, _batch: Optional[str] = None) -> str:
    batch_id = _require_batch(batch, _batch)
    return str(_tree_root() / f"{batch_id}.leaves")


def fn_nodes_sorted(batch: Optional[str] = None, _batch: Optional[str] = None) -> str:
    batch_id = _require_batch(batch, _batch)
    return str(_tree_root() / f"{batch_id}.nodes")


def fn_bubblejson(batch: Optional[str] = None, _batch: Optional[str] = None) -> str:
    batch_id = _require_batch(batch, _batch)
    return str(_bubble_root() / f"{batch_id}_k{KMER_LENGTH}_bubbles.json")


def fn_bubblelog(batch: Optional[str] = None, _batch: Optional[str] = None) -> str:
    batch_id = _require_batch(batch, _batch)
    return str(_bubble_root() / f"{batch_id}_k{KMER_LENGTH}_bubbles.log")


def fn_deletionbubbles(batch: Optional[str] = None, _batch: Optional[str] = None) -> str:
    batch_id = _require_batch(batch, _batch)
    return str(_bubble_root() / f"{batch_id}_k{KMER_LENGTH}_deletion-bubbles.json")


def fn_deletionbubbles_chkpoint(batch: Optional[str] = None, _batch: Optional[str] = None) -> str:
    batch_id = _require_batch(batch, _batch)
    return str(INTERMEDIATE_DIR / "deletion-bubbles" / batch_id)


def fn_pcdbg(batch: Optional[str] = None, _batch: Optional[str] = None) -> str:
    batch_id = _require_batch(batch, _batch)
    return str(OUTPUT_DIR / f"{batch_id}_k{KMER_LENGTH}_pcdbg.gfa")


#############################################
# Wildcard helpers and common file rewrites #
#############################################


def w_sample_source(wildcards):
    return BATCH_SAMPLE_MAP[wildcards["batch"]][wildcards["sample"]]


def generate_file_list(input_list_fn, output_list_fn, filename_function):
    input_path = Path(input_list_fn)
    output_path = Path(output_list_fn)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with input_path.open() as src, output_path.open("w") as dest:
        for raw in src:
            entry = raw.strip()
            if not entry:
                continue
            target = filename_function(entry)
            rel_path = os.path.relpath(target, output_path.parent)
            dest.write(rel_path + "\n")


def load_list(fn) -> Iterable[str]:
    path = Path(fn)
    if not path.exists():
        print(f"File not found {fn}, using empty list")
        return []
    with path.open() as handle:
        return [line.strip() for line in handle]
