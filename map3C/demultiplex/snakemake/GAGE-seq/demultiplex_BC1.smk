import pandas as pd
import os
from glob import glob

run_info = pd.read_csv("run_config.csv")
# Get barcode names
cell_ids = []
barcodes = config["demultiplex_protocols"]["GAGE-seq"]["BC1"]["barcodes"]
with open(barcodes) as f:
    for line in f:
        line = line.strip()
        if line[0] == ">":
            cell_ids.append(line[1:])

def get_R1_fastq(wildcards):
    plate = wildcards.plate
    fastq = run_info.loc[plate]["R1"]
    return fastq

def get_R2_fastq(wildcards):
    plate = wildcards.plate
    fastq = run_info.loc[plate]["R2"]
    return fastq

rule all:
    input:
        expand("{plate}-BC1-{cell_id}/{plate}-BC1-{cell_id}_indexed_R1.fastq.gz", 
               plate=run_info.index, cell_id=cell_ids),
        expand("{plate}-BC1-{cell_id}/{plate}-BC1-{cell_id}_indexed_R2.fastq.gz", 
               plate=run_info.index, cell_id=cell_ids),
        expand("{plate}-BC1_demultiplex_stats.txt", plate=run_info.index),
        expand("{plate}-BC1_unknown_barcode_R1.fastq.gz", plate=run_info.index),
        expand("{plate}-BC1_unknown_barcode_R2.fastq.gz", plate=run_info.index)

rule demultiplex:
    input:
        R1=get_R1_fastq,
        R2=get_R2_fastq
    output:
        R1=expand("{{plate}}-BC1-{cell_id}/{{plate}}-BC1-{cell_id}_indexed_R1.fastq.gz", cell_id=cell_ids),
        R2=expand("{{plate}}-BC1-{cell_id}/{{plate}}-BC1-{cell_id}_indexed_R2.fastq.gz", cell_id=cell_ids),
        stats="{plate}-BC1_demultiplex_stats.txt",
        R1_untrimmed="{plate}-BC1_unknown_barcode_R1.fastq.gz",
        R2_untrimmed="{plate}-BC1_unknown_barcode_R2.fastq.gz"
    threads:
        2
    params:
        cutadapt_out_R1=lambda wildcards: f"{wildcards.plate}-BC1-{{name}}/{wildcards.plate}-BC1-{{name}}_indexed_R1.fastq.gz",
        cutadapt_out_R2=lambda wildcards: f"{wildcards.plate}-BC1-{{name}}/{wildcards.plate}-BC1-{{name}}_indexed_R2.fastq.gz",
        extra=config["demultiplex_protocols"]["GAGE-seq"]["BC1"]["cutadapt_params"]
    conda:
        "map3C_preprocess_cutadapt"
    wildcard_constraints:
        plate="[A-Za-z0-9-]+(?=-BC1)"
    shell:
        """
        cutadapt -j {threads} -Z {params.extra} -u 23 --no-indels --action=none -g ^file:{barcodes} \
            -o {params.cutadapt_out_R2} \
            -p {params.cutadapt_out_R1} \
            --untrimmed-output {output.R2_untrimmed} \
            --untrimmed-paired-output {output.R1_untrimmed} \
            {input.R2} {input.R1} > {output.stats}
        """