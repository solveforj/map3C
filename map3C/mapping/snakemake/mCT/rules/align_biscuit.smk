
rule align_biscuit:
    input:
        get_trimmed_r1_fastq, 
        get_trimmed_r2_fastq
    output:
        temp("{id}_biscuit.bam")
    threads: 
        10
    params:
        reference_path=config["dna"]["align"]["biscuit"]["reference_path"],
        extra=config["dna"]["align"]["biscuit"]["biscuit_params"]
    conda:
        "map3C_utils"
    retries: 3
    shell:
        """
        biscuit align -@ {threads} {params.extra} {params.reference_path} {input} \
            | samtools sort -o {output} -O BAM 
        """

rule bsconv:
    input:
        rules.align_biscuit.output
    output:
        temp("{id}_biscuit_bsconv.bam")
    threads:
        1
    params:
        reference_path=config["dna"]["align"]["biscuit"]["reference_path"]
    conda:
        "map3C_utils"
    shell:
        """
        biscuit bsconv {params.reference_path} {input} {output} 
        """

rule name_sort_bam:
    input:
        rules.bsconv.output
    output:
        temp("{id}_biscuit_nsort.bam")
    threads:
        10
    conda:
        "map3C_utils"
    shell:
        """
        samtools sort -@ {threads} -n -o {output} {input}
        """

rule remove_rna_reads:
    input:
        rules.name_sort_bam.output
    output:
        bam=temp("{id}_biscuit_contam_filtered.bam"),
        stats=temp("{id}_biscuit_contam_stats.txt")
    conda:
        "map3C_tools"
    params:
        out_prefix=lambda wildcards: f"{wildcards.id}_biscuit",
        extra=config["dna"]["align"]["contamination_filter_params"]
    shell:
        """
        map3C contamination-filter \
            --bam {input} \
            --out-prefix {params.out_prefix} \
            --mate-annotation flag \
            {params.extra}
        """

rule dna_contam_stats:
    input:
        stats=rules.remove_rna_reads.output.stats
    output:
        stats=temp("{id}_biscuit_mod_contam_stats.txt")
    threads:
        1
    run:
        with open(input["stats"]) as f:
            all_cols = ["biscuit_" + i for i in f.readline().strip().split()]
            all_vals = f.readline().strip().split()
        with open(output["stats"], "w") as f:
            f.write("\t".join(all_cols) + "\n")
            f.write("\t".join(all_vals) + "\n")
