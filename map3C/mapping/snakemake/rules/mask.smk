
rule mask:
    input:
        get_namesorted_analysis_bam
    output:
        (
            "{id}_map3C_masked.bam"
            if last_bam_step == "mask"
            else temp("{id}_map3C_masked.bam")
        ),
    params:
        out_prefix=lambda wildcards: f"{wildcards.id}",
        mate_annotation=('--mate-annotation qname ' 
                                if trim_output == "separate" and not joint_alignments 
                                else '--mate-annotation flag '),
        extra=config["read_analysis"]["mask"]["mask_params"]
    conda:
        "map3C_tools"
    threads:
        1
    shell:
        'map3C mask-overlaps '
        '--bam {input} '
        '--out-prefix {params.out_prefix} '
        '{params.mate_annotation} '
        '{params.extra} '

def get_namesorted_analysis_bam(wildcards):
    return f"{wildcards.id}_map3C_masked.bam"
