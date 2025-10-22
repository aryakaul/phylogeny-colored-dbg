##
## Tree inference
##
import os.path
from pathlib import Path


SCRIPTS_DIR = Path(workflow.basedir) / "scripts"


rule tree_postprocessing:
    """
    Get a cleaned tree and all auxiliary files
    """
    output:
        nwk=fn_tree_sorted(_batch="{batch}"),
        leaves=fn_leaves_sorted(_batch="{batch}"),
        nodes=fn_nodes_sorted(_batch="{batch}"),
    input:
        nwk=fn_tree_dirty(_batch="{batch}"),
    params:
        script=str(SCRIPTS_DIR / "postprocess_tree.py"),
    conda:
        "../envs/ete4.yml"
    shell:
        """
        {params.script} \\
            --standardize \\
            --midpoint-outgroup \\
            --name-internals \\
            --ladderize \\
            -l {output.leaves} \\
            -n {output.nodes} \\
            {input.nwk} {output.nwk}
        """


# rule cp_final_tree_for_post_output:
#    output:
#        nw=fn_post_output_tree(_batch="{batch}"),
#    input:
#        nw=fn_tree_sorted(_batch="{batch}"),
#    shell:
#        """
#        cp "{input.nw}" "{output.nw}"
#        """


rule symlink_nwk_tree:
    """
    Symlink a phylogenetic tree if possible (nwk)
    """
    output:
        nwk=fn_tree_dirty(_batch="{batch}"),
    input:
        nwk=dir_input() / "{batch}.nwk",
    params:
        relative_path=lambda wildcards, input, output: os.path.relpath(
            input.nwk, start=os.path.dirname(output.nwk)
        ),
    shell:
        """
        ln -sf {params.relative_path} {output.nwk}
        """


if not config["trees_required"]:

    ruleorder: symlink_nwk_tree > tree_newick_mashtree

    rule tree_newick_mashtree:
        """
        Infer a phylogenetic tree from the assemblies belonging to a given batch
        """
        output:
            nwk=fn_tree_dirty(_batch="{batch}"),
        input:
            w_batch_asms,
        threads: config["mashtree_threads"]
        params:
            k=config["mashtree_kmer_length"],
            s=config["mashtree_sketch_size"],
            t=min(int(config["mashtree_threads"]), workflow.cores),  # ensure that the number of cores for MashTree doesn't go too low
        conda:
            "../envs/mashtree.yaml"
        shell:
            """
            mashtree \\
                --numcpus {params.t} \\
                --kmerlength {params.k} \\
                --sketch-size {params.s} \\
                --seed 42  \\
                {input} \\
                | tee {output.nwk}
            """
