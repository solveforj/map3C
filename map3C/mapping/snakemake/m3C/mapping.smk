import pandas as pd

run_info = pd.read_csv("run_config.csv")
    
mode = config["general"]["mode"]

sort_done = config["contacts"]["sort"]["sort_protocol"] != "none"
dedup_done = config["contacts"]["dedup"]["dedup_protocol"] != "none"
lowcov_done = config["contacts"]["lowcov"]["lowcov_protocol"] != "none"
call_done = config["contacts"]["call"]["call_protocol"] != "none"
mask_done = config["read_analysis"]["mask"]["mask_protocol"] != "none"
coord_sort_bam_done = config["read_analysis"]["coord_sort_bam"]["coord_sort_bam_protocol"] != "none"
allc_done = config["read_analysis"]["allc"]["allc_protocol"] != "none"

bam_generated = "--no-output-bam" not in config["contacts"]["call"]["call_params"]
keep_highcov = config["contacts"]["lowcov"]["keep_highcov"]

if call_done:
    last_contacts_step = "call"
    filter_suffix = "map3C"    
if sort_done:
    last_contacts_step = "sort"
    filter_suffix = "map3C.srt"
if dedup_done:
    last_contacts_step = "dedup"
    filter_suffix = "map3C.srt.dedup"
if lowcov_done:
    last_contacts_step = "lowcov"
    filter_suffix = "map3C.srt.dedup.lcov"
    highcov = "{id}_map3C.srt.dedup.hcov.pairs.gz"

pairs = f"{{id}}_{filter_suffix}.pairs.gz"

if call_done:
    last_bam_step = "call"
    bam_suffix = "map3C"
if mask_done:
    last_bam_step = "mask"
    bam_suffix += "_masked"
if coord_sort_bam_done:
    last_bam_step = "sort"
    bam_suffix += "_sorted"

bam = f"{{id}}_{bam_suffix}.bam"

if mode == "bsdna":
    
    rule all:
        input:
            # QC stats
            expand("{id}_qc_stats.txt", id=run_info.index),
            # Alignments
            (expand(bam, id=run_info.index)
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
            (expand(bam, id=run_info.index)
             if bam_generated
             else []),
            # Contacts
            expand(pairs, id=run_info.index),
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
