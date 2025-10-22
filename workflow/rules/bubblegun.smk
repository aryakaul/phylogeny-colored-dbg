rule run_bubblegun:
    output:
        bubblejson=fn_bubblejson(_batch="{batch}"),
        bubblelog=fn_bubblelog(_batch="{batch}"),
    input:
        gfa=fn_cuttlefish_out(_batch="{batch}", _ext="gfa1"),
    conda:
        "../envs/bubblegun.yml"
    shell:
        """
        mkdir -p $(dirname {output.bubblejson})
        BubbleGun \\
            -g {input.gfa} \\
            --log_file {output.bubblelog} \\
            bchains \\
            --bubble_json {output.bubblejson}
        """
