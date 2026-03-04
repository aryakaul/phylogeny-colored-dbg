"""
Evolutionary stratification analysis.

Decomposes pangenome graph complexity by phylogenetic depth, revealing
when in evolutionary history different layers of variation were generated.
"""

from pathlib import Path


SCRIPTS_DIR = Path(workflow.basedir) / "scripts"


# ---------------------------------------------------------------------------
# Output path helpers
# ---------------------------------------------------------------------------

def fn_evolutionary_stratification(batch: str | None = None, _batch: str | None = None) -> str:
    """Output path for evolutionary stratification TSV."""
    batch_id = batch if batch is not None else _batch
    if batch_id is None:
        raise ValueError("Batch identifier is required.")
    label = parsimony_label()
    return str(
        dir_output() / "stratification" / f"{batch_id}_k{KMER_LENGTH}_{label}_evolutionary_stratification.tsv"
    )


def fn_evolutionary_stratification_plot(batch: str | None = None, _batch: str | None = None) -> str:
    """Output path for evolutionary stratification SVG plot."""
    batch_id = batch if batch is not None else _batch
    if batch_id is None:
        raise ValueError("Batch identifier is required.")
    label = parsimony_label()
    return str(
        dir_output() / "stratification" / f"{batch_id}_k{KMER_LENGTH}_{label}_evolutionary_stratification.svg"
    )


# ---------------------------------------------------------------------------
# Rule: Evolutionary stratification
# ---------------------------------------------------------------------------

rule evolutionary_stratification:
    """
    Decompose pangenome graph complexity by phylogenetic depth.
    Sweeps a depth threshold from root to tips, measuring structural
    properties (components, branching, sequence content) at each level.
    """
    input:
        sqldb=fn_sqldb(_batch="{batch}"),
        tree=fn_tree_sorted(_batch="{batch}"),
    output:
        tsv=fn_evolutionary_stratification(_batch="{batch}"),
        plot=fn_evolutionary_stratification_plot(_batch="{batch}"),
    params:
        script=str(SCRIPTS_DIR / "evolutionary_stratification"),
    conda:
        "../envs/ete4.yml"
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        {params.script} \\
            --db {input.sqldb} \\
            --tree {input.tree} \\
            --output {output.tsv} \\
            --plot {output.plot} \\
            -v
        """
