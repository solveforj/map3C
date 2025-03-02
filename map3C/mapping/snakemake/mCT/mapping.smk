import pandas as pd

run_info = pd.read_csv("run_config.csv")
    
mode = config["general"]["mode"]

dna_done = config["dna"]["run"] == True
rna_done = config["rna"]["run"] == True

rule all:
    input:
        # QC stats
        expand("{id}_qc_stats.txt", id=run_info.index),
        # DNA alignments
        (expand("{id}_biscuit_csort.bam", id=run_info.index)
        if dna_done
        else []),
        # Methylation
        (expand("{id}.allc.tsv.gz.tbi", id=run_info.index)
        if dna_done
        else []),
        (expand("{id}.allc.tsv.gz.count.csv", id=run_info.index)
        if dna_done
        else []),
        (expand("{id}.allc.tsv.gz", id=run_info.index)
        if dna_done
        else []),
        # Pairs
        (expand("{id}_biscuit_map3C.pairs.gz", id=run_info.index)
        if dna_done
        else []),
        # RNA alignments
        (expand("{id}_STAR_PE_dedup.bam", id=run_info.index)
        if rna_done
        else []),
        (expand("{id}_STAR_SE1_dedup.bam", id=run_info.index)
        if rna_done
        else []),
        (expand("{id}_STAR_SE2_dedup.bam", id=run_info.index)
        if rna_done
        else []),
        # RNA exon counts
        (expand("{id}_STAR_exon_counts.txt", id=run_info.index)
        if rna_done
        else []),
        # RNA gene counts
        (expand("{id}_STAR_gene_counts.txt", id=run_info.index)
        if rna_done
        else []),

include: "rules/preprocess.smk"
include: "rules/align.smk"
include: "rules/postprocess.smk"

# STATS: 
include: "rules/stats.smk"
