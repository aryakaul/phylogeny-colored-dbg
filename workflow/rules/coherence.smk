"""
Evolutionary persistence analysis.

Measures the spatial extent over which evolutionary history is conserved
along graph walks via strict and containment coherence decay curves.
"""

from pathlib import Path


SCRIPTS_DIR = Path(workflow.basedir) / "scripts"


# ---------------------------------------------------------------------------
# Output path helpers
# ---------------------------------------------------------------------------


def fn_evolutionary_persistence(
    batch: str | None = None, _batch: str | None = None
) -> str:
    """Output path for evolutionary persistence TSV."""
    batch_id = batch if batch is not None else _batch
    if batch_id is None:
        raise ValueError("Batch identifier is required.")
    label = parsimony_label()
    return str(
        dir_output()
        / "persistence"
        / f"{batch_id}_k{KMER_LENGTH}_{label}_evolutionary_persistence.tsv"
    )


def fn_evolutionary_persistence_plot(
    batch: str | None = None, _batch: str | None = None
) -> str:
    """Output path for evolutionary persistence SVG plot."""
    batch_id = batch if batch is not None else _batch
    if batch_id is None:
        raise ValueError("Batch identifier is required.")
    label = parsimony_label()
    return str(
        dir_output()
        / "persistence"
        / f"{batch_id}_k{KMER_LENGTH}_{label}_evolutionary_persistence.svg"
    )


def fn_evolutionary_persistence_per_seed(
    batch: str | None = None, _batch: str | None = None
) -> str:
    """Output path for per-seed persistence TSV."""
    batch_id = batch if batch is not None else _batch
    if batch_id is None:
        raise ValueError("Batch identifier is required.")
    label = parsimony_label()
    return str(
        dir_output()
        / "persistence"
        / f"{batch_id}_k{KMER_LENGTH}_{label}_per_seed_persistence.tsv"
    )


# ---------------------------------------------------------------------------
# Rule: Evolutionary persistence
# ---------------------------------------------------------------------------


rule evolutionary_persistence:
    """
    Measure how rapidly phylogenetic assignment changes with graph distance.
    Produces decay curves for strict and containment coherence vs. hop count.
    """
    input:
        sqldb=fn_sqldb(_batch="{batch}"),
        tree=fn_tree_sorted(_batch="{batch}"),
    output:
        tsv=fn_evolutionary_persistence(_batch="{batch}"),
        plot=fn_evolutionary_persistence_plot(_batch="{batch}"),
        per_seed=fn_evolutionary_persistence_per_seed(_batch="{batch}"),
    params:
        script=str(SCRIPTS_DIR / "evolutionary_persistence"),
        max_hops=20,
        samples=5000,
    threads: MAX_THREADS
    conda:
        "../envs/ete4.yml"
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        {params.script} \\
            --db {input.sqldb} \\
            --tree {input.tree} \\
            --max-hops {params.max_hops} \\
            --samples {params.samples} \\
            --jobs {threads} \\
            --output {output.tsv} \\
            --plot {output.plot} \\
            --per-seed-output {output.per_seed} \\
            -v
        """
