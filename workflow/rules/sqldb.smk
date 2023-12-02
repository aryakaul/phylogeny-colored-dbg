rule build_sqldb:
    output:
        sqldb=fn_sqldb(_batch="{batch}"),
    input:
        gfa=fn_cuttlefish_out(_batch="{batch}", _ext="gfa1"),
        unitigs_to_cuts=fn_unitig2cuts(_batch="{batch}"),
    conda:
        "../envs/pandas.yml"
    params:
        script=snakemake.workflow.srcdir("../scripts/gfa_to_sql"),
    shell:
        """
        mkdir -p $(dirname {output})
        {params.script} \\
            -g {input.gfa} \\
            -s {input.unitigs_to_cuts} \\
            -o {output} \\
            -vv
        """
