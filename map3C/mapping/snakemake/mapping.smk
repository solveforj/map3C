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
generate_sr = "--enzymeless-split-read-pairs" in config["contacts"]["filter"]["filter_params"]
generate_enzyme = "--enzyme-pairs" in config["contacts"]["filter"]["filter_params"]
generate_enzymeless = "--enzymeless-pairs" in config["contacts"]["filter"]["filter_params"]

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
    
if filter_done:
    last_contacts_step = "filter"
    split_reads = f"{{id}}_{filter_suffix}.flt.enzymeless_split_reads.pairs.gz"
    enzyme = f"{{id}}_{filter_suffix}.flt.enzyme.pairs.gz"
    enzymeless = f"{{id}}_{filter_suffix}.flt.enzymeless.pairs.gz"
    pairs = f"{{id}}_{filter_suffix}.flt.pairs.gz"

if mode == "bsdna":
    
    rule all:
        input:
            # QC stats
            expand("{id}_qc_stats.txt", id=run_info.index),
            # Alignments
            (expand("{id}_map3C.bam", id=run_info.index)
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
            # Enzyme
            (expand(enzyme, id=run_info.index)
             if generate_enzyme and filter_done
             else []),
            # Enzymeless
            (expand(enzymeless, id=run_info.index)
             if generate_enzymeless and filter_done
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
            (expand("{id}_map3C.bam", id=run_info.index)
             if bam_generated
             else []),
            # Contacts
            expand(pairs, id=run_info.index),
            # Split reads
            (expand(split_reads, id=run_info.index)
             if generate_sr and filter_done
             else []),
            # Enzyme
            (expand(enzyme, id=run_info.index)
             if generate_enzyme and filter_done
             else []),
            # Enzymeless
            (expand(enzymeless, id=run_info.index)
             if generate_enzymeless and filter_done
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
