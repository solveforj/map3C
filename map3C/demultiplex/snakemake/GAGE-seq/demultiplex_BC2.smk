import pandas as pd
import os
from glob import glob

run_info = pd.read_csv("run_config.csv")
# Get barcode names
cell_ids = []
barcodes = config["demultiplex_protocols"]["GAGE-seq"]["BC2"]["barcodes"]
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
        expand("{plate}-BC2-{cell_id}/{plate}-BC2-{cell_id}_R1.fastq.gz", 
               plate=run_info.index, cell_id=cell_ids),
        expand("{plate}-BC2-{cell_id}/{plate}-BC2-{cell_id}_R2.fastq.gz", 
               plate=run_info.index, cell_id=cell_ids),
        expand("{plate}-BC2_demultiplex_stats.txt", plate=run_info.index),
        expand("{plate}-BC2_unknown_barcode_R1.fastq.gz", plate=run_info.index),
        expand("{plate}-BC2_unknown_barcode_R2.fastq.gz", plate=run_info.index)

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
        R1=expand("{{plate}}-BC2-{cell_id}/{{plate}}-BC2-{cell_id}_R1.fastq.gz", cell_id=cell_ids),
        R2=expand("{{plate}}-BC2-{cell_id}/{{plate}}-BC2-{cell_id}_R2.fastq.gz", cell_id=cell_ids),
        stats="{plate}-BC2_demultiplex_stats.txt",
        R1_untrimmed="{plate}-BC2_unknown_barcode_R1.fastq.gz",
        R2_untrimmed="{plate}-BC2_unknown_barcode_R2.fastq.gz"
    threads:
        2
    params:
        cutadapt_out_R1=lambda wildcards: f"{wildcards.plate}-BC2-{{name}}/{wildcards.plate}-BC2-{{name}}_R1.fastq.gz",
        cutadapt_out_R2=lambda wildcards: f"{wildcards.plate}-BC2-{{name}}/{wildcards.plate}-BC2-{{name}}_R2.fastq.gz",
        extra=config["demultiplex_protocols"]["GAGE-seq"]["BC2"]["cutadapt_params"]
    conda:
        "map3C_preprocess_cutadapt"
    wildcard_constraints:
        plate="[A-Za-z0-9-]+(?=-BC2)"
    shell:
        """
        cutadapt -j {threads} -Z {params.extra} --no-indels --action=none -g ^file:{barcodes} \
            -o {params.cutadapt_out_R2} \
            -p {params.cutadapt_out_R1} \
            --untrimmed-output {output.R2_untrimmed} \
            --untrimmed-paired-output {output.R1_untrimmed} \
            {input.R2} {input.R1} > {output.stats}
        """