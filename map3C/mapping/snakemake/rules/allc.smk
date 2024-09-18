
if mask_protocol == "none":

    rule coord_sort_trimmed:
        input:
            rules.generate_contacts.output.bam
        output:
            temp("{id}_trimmed_sorted.bam")
        threads: 
            10
        conda:
            "map3C_utils"
        shell:
            """
            samtools sort -@ {threads} -o {output} {input}
            """

    def get_coordsorted_analysis_bam(wildcards):
        return f"{wildcards.id}_trimmed_sorted.bam"

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
