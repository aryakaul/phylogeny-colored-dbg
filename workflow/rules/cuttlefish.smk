rule run_cuttlefish:
    output:
        gfa=fn_cuttlefish_out(_batch="{batch}", _ext="gfa1"),
        json=fn_cuttlefish_out(_batch="{batch}", _ext="json"),
    input:
        fof=f"{dir_input()}/{batch}.txt",
    conda:
        "../envs/cuttlefish.yml"
    params:
        k=config["kmer_length"],
    shell:
        """
        ulimit -n 2048
        OUTPUTDIR=$(dirname {output.gfa})
        mkdir -p $OUTPUTDIR
        cuttlefish build \
                -l {input} \
                -k {params.k} \
                -f 1 \
                -o $(dirname {output.gfa})/$(basename {output.gfa} .gfa1)
        """


rule unitig_sample_matrix_gfa1:
    output:
        binarymatrix=fn_colormtx(_batch="{batch}"),
    input:
        gfa=fn_cuttlefish_out(_batch="{batch}", _ext="gfa1"),
        fof=f"{dir_input()}/{batch}.txt",
    params:
        script=snakemake.workflow.srcdir(
            "../scripts/generate_unitig_colormatrix_cuttlefish_gfa1"
        ),
    conda:
        "../envs/pandas.yml"
    shell:
        """
        {params.script} \\
            {input.gfa} \\
            {input.fof} \\
            {output}
        """


rule mergeredundantcolors:
    output:
        uqcolors=fn_uqcolors(_batch="{batch}"),
        redundantcolors=fn_redundantcolors(_batch="{batch}"),
    input:
        binarymatrix=fn_colormtx(_batch="{batch}"),
    params:
        script=snakemake.workflow.srcdir("../scripts/uqcolorvectors_from_mtx"),
    conda:
        "../envs/pandas.yml"
    shell:
        """
        {params.script} \\
            {input.binarymatrix} \\
            {output.uqcolors} \\
            {output.redundantcolors}
        """
