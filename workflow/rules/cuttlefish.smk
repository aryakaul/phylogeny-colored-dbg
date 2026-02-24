from pathlib import Path


SCRIPTS_DIR = Path(workflow.basedir) / "scripts"


rule run_cuttlefish:
    output:
        gfa=fn_cuttlefish_out(_batch="{batch}", _ext="gfa1"),
        json=fn_cuttlefish_out(_batch="{batch}", _ext="json"),
    input:
        fof=f"{dir_input()}" + "/{batch}.txt",
    conda:
        "../envs/cuttlefish.yml"
    params:
        k=config["kmer_length"],
        temp=f"{config.get('tmp_dir', 'tmp')}" + "/{batch}_cuttlefish",
    threads: MAX_THREADS
    shell:
        """
        ulimit -n 2048
        OUTPUTDIR=$(dirname {output.gfa})
        mkdir -p {params.temp}
        mkdir -p $OUTPUTDIR
        cuttlefish build \
                -l {input} \
                -k {params.k} \
                -t {threads} \
                -f 1 \
                --unrestrict-memory \
                -w {params.temp} \
                -o $(dirname {output.gfa})/$(basename {output.gfa} .gfa1)
        rm -r {params.temp}
        """


rule unitig_sample_matrix_gfa1:
    output:
        binarymatrix=fn_colormtx(_batch="{batch}"),
    input:
        gfa=fn_cuttlefish_out(_batch="{batch}", _ext="gfa1"),
        fof=f"{dir_input()}" + "/{batch}.txt",
    params:
        script=str(SCRIPTS_DIR / "generate_unitig_colormatrix_cuttlefish_gfa1"),
    conda:
        "../envs/pandas.yml"
    shell:
        """
        {params.script} -v \\
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
        script=str(SCRIPTS_DIR / "uqcolorvectors_from_mtx"),
    conda:
        "../envs/pandas.yml"
    shell:
        """
        {params.script} \\
            {input.binarymatrix} \\
            {output.uqcolors} \\
            {output.redundantcolors}
        """
