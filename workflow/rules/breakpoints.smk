"""
Concordance breakpoint detection rules.

Detects edges in the pcDBG where adjacent unitigs have discordant
phylogenetic assignments, indicating boundaries between differently-inherited
genomic regions.
"""

from pathlib import Path


SCRIPTS_DIR = Path(workflow.basedir) / "scripts"


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
        min_distance=1,  # Minimum tree distance to report
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

        # Tree distance distribution
        if len(df) > 0:
            dist_counts = df["tree_distance"].value_counts().sort_index()
            for dist, count in dist_counts.items():
                summary[f"distance_{dist}_count"] = count

        # Write summary
        with open(output.summary, "w") as f:
            for key, value in summary.items():
                f.write(f"{key}\t{value}\n")
