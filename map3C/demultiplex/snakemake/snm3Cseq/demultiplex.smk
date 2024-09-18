import pandas as pd
import os
from glob import glob

run_info = pd.read_csv("run_config.csv")
# Get barcode names
cell_ids = []
barcodes = config["demultiplex_protocols"]["snm3Cseq"]["barcodes"]
with open(barcodes) as f:
    for line in f:
        line = line.strip()
        if line[0] == ">":
            cell_ids.append(line[1:])

def get_R1_fastq(wildcards):
    plate = wildcards.plate
    fastq_dir = run_info.loc[plate]["fastq_dir"]
    fastq_files = sorted(glob(os.path.join(fastq_dir, f"*{plate}[-_]*R1*.fastq.gz")))
    return fastq_files

def get_R2_fastq(wildcards):
    plate = wildcards.plate
    fastq_dir = run_info.loc[plate]["fastq_dir"]
    fastq_files = sorted(glob(os.path.join(fastq_dir, f"*{plate}[-_]*R2*.fastq.gz")))
    return fastq_files

rule all:
    input:
        expand("{plate}_{cell_id}/{plate}_{cell_id}_indexed_R1.fastq.gz", 
               plate=run_info.index, cell_id=cell_ids),
        expand("{plate}_{cell_id}/{plate}_{cell_id}_indexed_R2.fastq.gz", 
               plate=run_info.index, cell_id=cell_ids),
        expand("{plate}_demultiplex_stats.txt", plate=run_info.index),
        expand("{plate}_unknown_barcode_R1.fastq.gz", plate=run_info.index),
        expand("{plate}_unknown_barcode_R2.fastq.gz", plate=run_info.index)

rule concatenate_R1:
    input:
        get_R1_fastq
    output:
        temp("{plate}_R1.fastq.gz")
    conda:
        "map3C_utils"
    wildcard_constraints:
        plate="[A-Za-z0-9-]+(?=_|\\.)"
    shell:
        'cat {input} > {output}'

rule concatenate_R2:
    input:
        get_R2_fastq
    output:
        temp("{plate}_R2.fastq.gz")
    conda:
        "map3C_utils"
    wildcard_constraints:
        plate="[A-Za-z0-9-]+(?=_|\\.)"
    shell:
        'cat {input} > {output}'

rule demultiplex:
    input:
        R1=rules.concatenate_R1.output,
        R2=rules.concatenate_R2.output
    output:
        R1=expand("{{plate}}_{cell_id}/{{plate}}_{cell_id}_indexed_R1.fastq.gz", cell_id=cell_ids),
        R2=expand("{{plate}}_{cell_id}/{{plate}}_{cell_id}_indexed_R2.fastq.gz", cell_id=cell_ids),
        stats="{plate}_demultiplex_stats.txt",
        R1_untrimmed="{plate}_unknown_barcode_R1.fastq.gz",
        R2_untrimmed="{plate}_unknown_barcode_R2.fastq.gz"
    threads:
        2
    params:
        cutadapt_out_R1=lambda wildcards: f"{wildcards.plate}_{{name}}/{wildcards.plate}_{{name}}_indexed_R1.fastq.gz",
        cutadapt_out_R2=lambda wildcards: f"{wildcards.plate}_{{name}}/{wildcards.plate}_{{name}}_indexed_R2.fastq.gz",
        extra=config["demultiplex_protocols"]["snm3Cseq"]["cutadapt_params"]
    conda:
        "map3C_preprocess_cutadapt"
    wildcard_constraints:
        plate="[A-Za-z0-9-]+(?=_|\\.)"
    shell:
        """
        cutadapt -j {threads} -Z {params.extra} --no-indels --action=none -g ^file:{barcodes} \
            -o {params.cutadapt_out_R1} \
            -p {params.cutadapt_out_R2} \
            --untrimmed-output {output.R1_untrimmed} \
            --untrimmed-paired-output {output.R2_untrimmed} \
            {input.R1} {input.R2} > {output.stats}
        """