from pathlib import Path


SCRIPTS_DIR = Path(workflow.basedir) / "scripts"


rule add_colorinfo_to_gfa:
    input:
        sqldb=fn_sqldb(_batch="{batch}"),
        gfa=fn_cuttlefish_out(_batch="{batch}", _ext="gfa1"),
        unitigs_to_cuts=fn_unitig2cuts(_batch="{batch}"),
    output:
        pcdbg=fn_pcdbg(_batch="{batch}"),
    params:
        script=str(SCRIPTS_DIR / "add_coloredgfatag_cuttlefish-gfa1"),
        maxcolors=70,
    conda:
        "../envs/pandas.yml"
    shell:
        """
        mkdir -p $(dirname {output})
        {params.script} \\
            -g {input.gfa} \\
            -m {input.unitigs_to_cuts} \\
            -d {input.sqldb} \\
            -o {output} \\
            -vv \\
            -mc {params.maxcolors}
        """
