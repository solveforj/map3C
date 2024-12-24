
rule bam_to_allc:
    input:
        get_coordsorted_analysis_bam
    output:
        allc = "{id}.allc.tsv.gz",
        tbi = "{id}.allc.tsv.gz.tbi",
        stats = "{id}.allc.tsv.gz.count.csv",
        methylation_stats = temp("{id}_methylation_stats.txt")
    params:
        out_prefix=lambda wildcards: f"{wildcards.id}",
        reference_fasta=config["align"]["align_params"]["biscuit"]["reference_path"],
        extra=config["read_analysis"]["allc"]["allc_params"]
    conda:
        "map3C_tools"
    threads:
        1
    shell:
        'map3C bam-to-allc '
        '--bam-path {input} '
        '--reference-fasta {params.reference_fasta} '
        '--out-prefix {params.out_prefix} ' 
        '--save-count-df ' 
        '{params.extra} '
