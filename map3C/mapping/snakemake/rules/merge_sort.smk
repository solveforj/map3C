
rule merge:
    input:
        get_bams
    output:
        temp("{id}_merged.bam")
    threads: 
        10
    conda:
        "map3C_utils"
    shell:
        """
        samtools merge -@ {threads} -o {output} {input}
        """

rule sort:
    input:
        rules.merge.output
    output:
        temp("{id}_merged_sorted.bam")
    threads: 
        10
    conda:
        "map3C_utils"
    shell:
        """
        samtools sort -@ {threads} -n -o {output} {input}
        """

contamination_protocol = config["contamination"]["contamination_protocol"]

if mode == "bsdna":

    if contamination_protocol == "default": 

        rule contam_filter:
            input:
                bam = rules.sort.output
            output:
                bam=temp("{id}_contam_filtered.bam"),
                stats=temp("{id}_contam_stats.txt")      
            params:
                out_prefix=lambda wildcards: f"{wildcards.id}",
                mate_annotation=('--mate-annotation qname ' 
                                if trim_output == "separate" and not joint_alignments 
                                else '--mate-annotation flag '),
                extra = config["contamination"]["params"]
            conda:
                "map3C_tools"
            threads:
                1
            shell:
                'map3C contamination-filter '
                '{params.mate_annotation}  '
                '--bam {input.bam} '
                '--out-prefix {params.out_prefix} '
                '{params.extra} '
    
    
        def get_merged_bam(wildcards):
            return f"{wildcards.id}_contam_filtered.bam"

    else:

        def get_merged_bam(wildcards):
            return f"{wildcards.id}_merged_sorted.bam"

else:

    def get_merged_bam(wildcards):
        return f"{wildcards.id}_merged_sorted.bam"
