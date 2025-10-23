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

    ruleorder: symlink_nwk_tree > tree_newick_attotree

    rule tree_newick_attotree:
        """
        Infer a phylogenetic tree for the batch using attotree (distance + NJ).
        """
        output:
            nwk=fn_tree_dirty(_batch="{batch}"),
            distance=fn_tree_distance(_batch="{batch}"),
        input:
            manifest=dir_input() / "{batch}.txt",
        conda:
            "../envs/attotree.yml"
        shell:
            """
            mkdir -p $(dirname {output.nwk})
            mkdir -p $(dirname {output.distance})
            attotree \\
                -L {input.manifest} \\
                -o {output.nwk} \\
                -D \\
            > {output.distance}
            """
