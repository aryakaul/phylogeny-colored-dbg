"""
Phylogenetic coherence decay analysis.

Measures how rapidly phylogenetic assignment (cut-string) changes as
you walk outward from a unitig in the pcDBG. Produces a decay curve
and optional plot.
"""

from pathlib import Path


SCRIPTS_DIR = Path(workflow.basedir) / "scripts"


# ---------------------------------------------------------------------------
# Output path helpers
# ---------------------------------------------------------------------------

def fn_coherence_decay(batch: str | None = None, _batch: str | None = None) -> str:
    """Output path for phylogenetic coherence decay TSV."""
    batch_id = batch if batch is not None else _batch
    if batch_id is None:
        raise ValueError("Batch identifier is required.")
    label = parsimony_label()
    return str(
        dir_output() / "coherence" / f"{batch_id}_k{KMER_LENGTH}_{label}_coherence_decay.tsv"
    )


def fn_coherence_decay_plot(batch: str | None = None, _batch: str | None = None) -> str:
    """Output path for coherence decay PNG plot."""
    batch_id = batch if batch is not None else _batch
    if batch_id is None:
        raise ValueError("Batch identifier is required.")
    label = parsimony_label()
    return str(
        dir_output() / "coherence" / f"{batch_id}_k{KMER_LENGTH}_{label}_coherence_decay.png"
    )


# ---------------------------------------------------------------------------
# Rule: Phylogenetic coherence decay
# ---------------------------------------------------------------------------

rule phylogenetic_coherence_decay:
    """
    Measure how rapidly phylogenetic assignment changes with graph distance.
    Produces a decay curve: fraction of same-colored neighbours vs. hop count.
    """
    input:
        sqldb=fn_sqldb(_batch="{batch}"),
    output:
        tsv=fn_coherence_decay(_batch="{batch}"),
        plot=fn_coherence_decay_plot(_batch="{batch}"),
    params:
        script=str(SCRIPTS_DIR / "phylogenetic_coherence_decay"),
        max_hops=20,
        samples=10000,
    conda:
        "../envs/ete4.yml"
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        {params.script} \\
            --db {input.sqldb} \\
            --max-hops {params.max_hops} \\
            --samples {params.samples} \\
            --output {output.tsv} \\
            --plot {output.plot} \\
            -v
        """
