import pandas as pd

run_info = pd.read_csv("run_config.csv")
    
mode = config["general"]["mode"]

sort_done = config["contacts"]["sort"]["sort_protocol"] != "none"
dedup_done = config["contacts"]["dedup"]["dedup_protocol"] != "none"
lowcov_done = config["contacts"]["lowcov"]["lowcov_protocol"] != "none"
call_done = config["contacts"]["call"]["call_protocol"] != "none"
filter_done = config["contacts"]["filter"]["filter_protocol"] != "none"
allc_done = config["read_analysis"]["allc"]["allc_protocol"]

bam_generated = "--no-output-bam" not in config["contacts"]["call"]["call_params"]
keep_highcov = config["contacts"]["lowcov"]["keep_highcov"]
generate_sr = "--split-reads" in config["contacts"]["filter"]["filter_params"]

if call_done:
    last_contacts_step = "call"
    filter_suffix = "all"    
if sort_done:
    last_contacts_step = "sort"
    filter_suffix = "all.srt"
if dedup_done:
    last_contacts_step = "dedup"
    filter_suffix = "all.srt.dedup"
if lowcov_done:
    last_contacts_step = "lowcov"
    filter_suffix = "all.srt.dedup.lcov"
    highcov = "{id}_artefacts.srt.dedup.hcov.pairs.gz"

pairs = f"{{id}}_{filter_suffix}.pairs.gz"
    
if filter_done:
    split_reads = pairs.replace(".pairs.gz", ".split_reads.pairs.gz")
    pairs = pairs.replace(".pairs.gz", ".flt.pairs.gz")

if mode == "bsdna":
    
    rule all:
        input:
            # QC stats
            expand("{id}_qc_stats.txt", id=run_info.index),
            # Alignments
            (expand("{id}_trimmed.bam", id=run_info.index)
             if bam_generated
             else []),
            # Methylation
            (expand("{id}.allc.tsv.gz.tbi", id=run_info.index)
             if allc_done
             else []),
            (expand("{id}.allc.tsv.gz.count.csv", id=run_info.index)
             if allc_done
             else []),
            (expand("{id}.allc.tsv.gz", id=run_info.index)
             if allc_done
             else []),
            # Contacts
            expand(pairs, id=run_info.index),
            # Split reads
            (expand(split_reads, id=run_info.index)
             if generate_sr and filter_done
             else []),
            # Highcov artefacts
            (expand(highcov, id=run_info.index)
             if keep_highcov and lowcov_done
             else [])
            
if mode == "dna":
    
    rule all:
        input:
            # QC stats
            expand("{id}_qc_stats.txt", id=run_info.index),
            # Alignments
            (expand("{id}_trimmed.bam", id=run_info.index)
             if bam_generated
             else []),
            # Contacts
            expand(pairs, id=run_info.index),
            # Split reads
            (expand(split_reads, id=run_info.index)
             if generate_sr and filter_done
             else []),
            # Highcov artefacts
            (expand(highcov, id=run_info.index)
             if keep_highcov and lowcov_done
             else [])

include: "rules/preprocess.smk"
include: "rules/reformat.smk"
include: "rules/align.smk"
include: "rules/merge_sort.smk"
include: "rules/mkdup.smk"
include: "rules/contacts.smk"
include: "rules/read_analysis.smk"
include: "rules/stats.smk"
