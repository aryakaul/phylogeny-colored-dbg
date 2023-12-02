rule find_colored_deletions:
    output:
        deletionbubbles=fn_deletionbubbles(_batch="{batch}"),
    input:
        bubblejson=fn_bubblejson(_batch="{batch}"),
        sql=fn_sqldb(_batch="{batch}"),
    params:
        script=snakemake.workflow.srcdir("../scripts/find_colored_deletions"),
        predeletionsize=1000,
        postdeletionsize=100,
        maxpaths=100000,
        timeout=600,
    conda:
        "../envs/pandas.yml"
    shell:
        """
        {params.script} \\
            -b {input.bubblejson} \\
            --db {input.sql} \\
            -o {output} \\
            -vv \\
            --maxpaths {params.maxpaths} \\
            -to {params.timeout} \\
            --predeletionsize {params.predeletionsize} \\
            --postdeletionsize {params.postdeletionsize} \\
        """


checkpoint find_paths_for_deletions:
    output:
        directory(fn_deletionbubbles_chkpoint(_batch="{batch}")),
    input:
        sql=fn_sqldb(_batch="{batch}"),
        deletionbubbles=fn_deletionbubbles(_batch="{batch}"),
    params:
        script=snakemake.workflow.srcdir("../scripts/paths_to_seqs"),
        kmer_length=config["kmer_length"],
    conda:
        "../envs/pandas.yml"
    shell:
        """
        {params.script} \\
            -d {input.sql} \\
            -k {params.kmer_length} \\
            -b {input.deletionbubbles} \\
            -vv \\
            -o {output}
        """
