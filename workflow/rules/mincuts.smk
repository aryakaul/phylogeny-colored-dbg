rule minimal_cuts_percolorset:
    output:
        cuts=fn_minimalcuts(_batch="{batch}"),
    input:
        tree=fn_tree_sorted(_batch="{batch}"),
        uqcolors=fn_uqcolors(_batch="{batch}"),
    params:
        script=snakemake.workflow.srcdir("../scripts/minimal_cuts"),
        plotdir=fn_minimalcuts_plotdir(_batch="{batch}"),
        chunksize=5000,
    conda:
        "../envs/ete3.yml"
    shell:
        """
        mkdir -p $(dirname {output})
        {params.script} \\
            {input.tree} \\
            {input.uqcolors} \\
            -o {output} \\
            -v \\
            --all \\
            -cs {params.chunksize}
        """


rule decompose_redundant_colors:
    output:
        unitigs_to_cuts=fn_unitig2cuts(_batch="{batch}"),
    input:
        cuts=fn_minimalcuts(_batch="{batch}"),
        redundantcolors=fn_redundantcolors(_batch="{batch}"),
    params:
        script=snakemake.workflow.srcdir("../scripts/map_unitigs-to-cuts"),
    conda:
        "../envs/pandas.yml"
    shell:
        """
        {params.script} \\
            {input.redundantcolors} \\
            {input.cuts} \\
            {output}
        """
