
if mask_done:

    coord_sort_bam = "{id}_map3C_masked_sorted.bam"

    def get_coordsorted_analysis_bam(wildcards):
        return f"{wildcards.id}_map3C_masked_sorted.bam"

else:
    coord_sort_bam = "{id}_map3C_sorted.bam"

    def get_coordsorted_analysis_bam(wildcards):
        return f"{wildcards.id}_map3C_sorted.bam"

rule coord_sort:
    input:
        get_namesorted_analysis_bam
    output:
        coord_sort_bam
    threads: 
        10
    conda:
        "map3C_utils"
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        """


