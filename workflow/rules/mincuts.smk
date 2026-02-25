from pathlib import Path


SCRIPTS_DIR = Path(workflow.basedir) / "scripts"


_USE_FITCH_SCRIPT = parsimony_mode() == "fitch"


rule minimal_cuts_percolorset:
    output:
        cuts=fn_minimalcuts(_batch="{batch}"),
        colors=fn_minimalcuts_colors(_batch="{batch}") if PRODUCE_MINIMALCUTS_COLORS else [],
    input:
        tree=fn_tree_sorted(_batch="{batch}"),
        uqcolors=fn_uqcolors(_batch="{batch}"),
    params:
        script=str(SCRIPTS_DIR / ("minimal_cuts_fitch" if _USE_FITCH_SCRIPT else "minimal_cuts")),
        chunksize=5000,
        mode=parsimony_mode(),
        branch_penalty_flag=parsimony_branch_penalty_flag(),
        colors_flag=lambda wc: f"-c {fn_minimalcuts_colors(batch=wc.batch)}" if PRODUCE_MINIMALCUTS_COLORS else "",
        use_fitch=_USE_FITCH_SCRIPT,
    conda:
        "../envs/ete4.yml"
    threads: MAX_THREADS
    shell:
        """

        mkdir -p $(dirname {output.cuts})
        if [ "{params.use_fitch}" = "True" ]; then
            {params.script} \\
                {input.tree} \\
                {input.uqcolors} \\
                -o {output.cuts} \\
                {params.colors_flag} \\
                -v \\
                --all \\
                -cs {params.chunksize}
        else
            {params.script} \\
                {input.tree} \\
                {input.uqcolors} \\
                -o {output.cuts} \\
                {params.colors_flag} \\
                -v \\
                --all \\
                --mode {params.mode} \\
                {params.branch_penalty_flag} \\
                -cs {params.chunksize} \\
                -j {threads}
        fi
        """


rule decompose_redundant_colors:
    output:
        unitigs_to_cuts=fn_unitig2cuts(_batch="{batch}"),
    input:
        cuts=fn_minimalcuts(_batch="{batch}"),
        redundantcolors=fn_redundantcolors(_batch="{batch}"),
    params:
        script=str(SCRIPTS_DIR / "map_unitigs-to-cuts"),
    conda:
        "../envs/pandas.yml"
    shell:
        """
        {params.script} \\
            {input.redundantcolors} \\
            {input.cuts} \\
            {output}
        """
