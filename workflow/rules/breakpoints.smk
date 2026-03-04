"""
Concordance breakpoint detection.

Rules for:
  1. Binary breakpoint detection (original concordance_breakpoints)
  2. Stratified breakpoints with permutation null model
  3. Summary statistics
"""

from pathlib import Path


SCRIPTS_DIR = Path(workflow.basedir) / "scripts"


# ---------------------------------------------------------------------------
# Output path helpers
# ---------------------------------------------------------------------------

def fn_breakpoints(batch: str | None = None, _batch: str | None = None) -> str:
    """Output path for concordance breakpoints TSV."""
    batch_id = batch if batch is not None else _batch
    if batch_id is None:
        raise ValueError("Batch identifier is required.")
    label = parsimony_label()
    return str(
        dir_output() / "breakpoints" / f"{batch_id}_k{KMER_LENGTH}_{label}_breakpoints.tsv"
    )


def fn_breakpoint_summary(batch: str | None = None, _batch: str | None = None) -> str:
    """Output path for breakpoint summary statistics."""
    batch_id = batch if batch is not None else _batch
    if batch_id is None:
        raise ValueError("Batch identifier is required.")
    label = parsimony_label()
    return str(
        dir_output() / "breakpoints" / f"{batch_id}_k{KMER_LENGTH}_{label}_breakpoint_summary.tsv"
    )


def fn_breakpoints_stratified(batch: str | None = None, _batch: str | None = None) -> str:
    """Output path for stratified concordance breakpoints TSV."""
    batch_id = batch if batch is not None else _batch
    if batch_id is None:
        raise ValueError("Batch identifier is required.")
    label = parsimony_label()
    return str(
        dir_output() / "breakpoints" / f"{batch_id}_k{KMER_LENGTH}_{label}_breakpoints_stratified.tsv"
    )


def fn_breakpoints_stratified_plot(batch: str | None = None, _batch: str | None = None) -> str:
    """Output path for stratified concordance breakpoints SVG plot."""
    batch_id = batch if batch is not None else _batch
    if batch_id is None:
        raise ValueError("Batch identifier is required.")
    label = parsimony_label()
    return str(
        dir_output() / "breakpoints" / f"{batch_id}_k{KMER_LENGTH}_{label}_breakpoints_stratified.svg"
    )


# ---------------------------------------------------------------------------
# Rule 1: Original binary breakpoint detection
# ---------------------------------------------------------------------------

rule detect_concordance_breakpoints:
    """
    Detect edges where adjacent unitigs have discordant phylogenetic assignments.
    """
    input:
        sqldb=fn_sqldb(_batch="{batch}"),
        tree=fn_tree_sorted(_batch="{batch}"),
    output:
        breakpoints=fn_breakpoints(_batch="{batch}"),
    params:
        script=str(SCRIPTS_DIR / "concordance_breakpoints"),
        min_distance=1,
    conda:
        "../envs/ete4.yml"
    shell:
        """
        mkdir -p $(dirname {output.breakpoints})
        {params.script} \\
            -d {input.sqldb} \\
            -t {input.tree} \\
            -o {output.breakpoints} \\
            --min-distance {params.min_distance} \\
            -v
        """


rule summarize_breakpoints:
    """
    Compute summary statistics for detected breakpoints.
    """
    input:
        breakpoints=fn_breakpoints(_batch="{batch}"),
    output:
        summary=fn_breakpoint_summary(_batch="{batch}"),
    run:
        import pandas as pd

        df = pd.read_csv(input.breakpoints, sep="\t")

        summary = {
            "total_breakpoints": len(df),
            "unique_unitigs_a": df["unitig_a"].nunique(),
            "unique_unitigs_b": df["unitig_b"].nunique(),
            "mean_tree_distance": df["tree_distance"].mean() if len(df) > 0 else 0,
            "max_tree_distance": df["tree_distance"].max() if len(df) > 0 else 0,
            "unique_assignments_a": df["assignment_a"].nunique() if len(df) > 0 else 0,
            "unique_assignments_b": df["assignment_b"].nunique() if len(df) > 0 else 0,
        }

        if len(df) > 0:
            dist_counts = df["tree_distance"].value_counts().sort_index()
            for dist, count in dist_counts.items():
                summary[f"distance_{dist}_count"] = count

        with open(output.summary, "w") as f:
            for key, value in summary.items():
                f.write(f"{key}\t{value}\n")


# ---------------------------------------------------------------------------
# Rule 2: Stratified breakpoints with permutation null model
# ---------------------------------------------------------------------------

rule stratified_concordance_breakpoints:
    """
    Continuous tree-distance metric per edge with permutation p-values.
    Replaces binary concordant/discordant with a distance + null model.
    """
    input:
        sqldb=fn_sqldb(_batch="{batch}"),
        tree=fn_tree_sorted(_batch="{batch}"),
    output:
        breakpoints=fn_breakpoints_stratified(_batch="{batch}"),
        plot=fn_breakpoints_stratified_plot(_batch="{batch}"),
    params:
        script=str(SCRIPTS_DIR / "concordance_breakpoints_stratified"),
        permutations=1000,
    threads: MAX_THREADS
    conda:
        "../envs/ete4.yml"
    shell:
        """
        mkdir -p $(dirname {output.breakpoints})
        {params.script} \\
            --db {input.sqldb} \\
            --tree {input.tree} \\
            --permutations {params.permutations} \\
            --jobs {threads} \\
            --output {output.breakpoints} \\
            --plot {output.plot} \\
            -v
        """
