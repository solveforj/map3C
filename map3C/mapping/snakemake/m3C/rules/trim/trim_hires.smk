
import subprocess
import shlex

rule get_dna_reads:
    input:
        unpack(get_fastq)
    output:
        dna_r1=temp("{id}_DNA_R1.fastq.gz"),
        dna_r2=temp("{id}_DNA_R2.fastq.gz"),
        rna_r1=temp("{id}_RNA_R1.fastq.gz"),
        rna_r2=temp("{id}_RNA_R2.fastq.gz"),
        stats=temp("{id}_cutadapt_stats.txt")
    threads:
        2
    params:
        extra=config["trim_methods"]["hires"]["cutadapt_params"]
    conda:
        "map3C_preprocess_hires"
    shell:
        """
        cutadapt --report=minimal {params.extra} -j {threads} \
            --untrimmed-output {output.dna_r1} --untrimmed-paired-output {output.dna_r2} \
            -o {output.rna_r1} -p {output.rna_r2} \
            {input.r1} {input.r2} > {output.stats}
        """

rule seqtk_stats:
    input:
        r1=rules.get_dna_reads.output.dna_r1
    output:
        stats=temp("{id}_seqtk_stats.txt")
    threads:
        2
    conda:
        "map3C_preprocess_hires"
    shell:
        """
        seqtk size {input.r1} > {output.stats}
        """

rule trim_stats:
    input:
        cutadapt_stats=rules.get_dna_reads.output.stats,
        seqtk_stats=rules.seqtk_stats.output.stats,
    output:
        stats = temp("{id}_trim_stats.txt")
    threads:
        1
    run:
        with open(output["stats"], "w") as f, \
            open(input["cutadapt_stats"]) as cutadapt, \
            open(input["seqtk_stats"]) as seqtk:
                
            input_data = cutadapt.readlines()[1].strip().split()
            in_reads = input_data[1]
            out_rna_reads = input_data[6]

            input_data = seqtk.readlines()[0].strip().split()
            out_dna_reads = input_data[0]
    
            f.write("\t".join(["pre_trim_pairs", "post_trim_rna_pairs", "post_trim_dna_pairs"]) + "\n")
            f.write("\t".join([in_reads, out_rna_reads, out_dna_reads]) + "\n")
    
trim_output = "separate"

def get_trimmed_r1_fastq(wildcards):
    return f"{wildcards.id}_DNA_R1.fastq.gz"

def get_trimmed_r2_fastq(wildcards):
    return f"{wildcards.id}_DNA_R2.fastq.gz"
