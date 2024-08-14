import pandas as pd

run_info = pd.read_csv("run_config.csv")
    
mode = config["general"]["mode"]

sort_done = config["contacts"]["sort"]["sort_protocol"] != "none"
dedup_done = config["contacts"]["dedup"]["dedup_protocol"] != "none"
lowcov_done = config["contacts"]["lowcov"]["lowcov_protocol"] != "none"
call_done = config["contacts"]["call"]["call_protocol"] != "none"

if call_done:
    last_contacts_step = "call"
    contacts = "{id}_contacts.pairs.gz"
    artefacts = "{id}_artefacts.pairs.gz"
if sort_done:
    last_contacts_step = "sort"
    contacts = "{id}_contacts.sorted.pairs.gz"
    artefacts = "{id}_artefacts.sorted.pairs.gz"
if dedup_done:
    last_contacts_step = "dedup"
    contacts = "{id}_contacts.dedup.pairs.gz"
    artefacts = "{id}_artefacts.dedup.pairs.gz"
if lowcov_done:
    last_contacts_step = "lowcov"
    contacts = "{id}_contacts.dedup.lowcov.pairs.gz"
    artefacts = "{id}_artefacts.dedup.pairs.gz"

if mode == "bsdna":
    
    rule all:
        input:
            # QC stats
            expand("{id}_qc_stats.txt", id=run_info.index),
            # Alignments
            expand("{id}_trimmed.bam", id=run_info.index),
            # Methylation
            (expand("{id}.allc.tsv.gz.tbi", id=run_info.index)
             if config["read_analysis"]["allc"]["allc_protocol"] != "none"
             else []),
            (expand("{id}.allc.tsv.gz.count.csv", id=run_info.index)
             if config["read_analysis"]["allc"]["allc_protocol"] != "none"
             else []),
            (expand("{id}.allc.tsv.gz", id=run_info.index)
             if config["read_analysis"]["allc"]["allc_protocol"] != "none"
             else []),
            # Contacts
            expand(contacts, id=run_info.index),
            # Artefacts
            expand(artefacts, id=run_info.index)

if mode == "dna":
    
    rule all:
        input:
            # QC stats
            expand("{id}_qc_stats.txt", id=run_info.index),
            # Alignments
            expand("{id}_trimmed.bam", id=run_info.index),
            # Contacts
            expand(contacts, id=run_info.index),
            # Artefacts
            expand(artefacts, id=run_info.index)


include: "rules/preprocess.smk"
include: "rules/reformat.smk"
include: "rules/align.smk"
include: "rules/merge_sort.smk"
include: "rules/mkdup.smk"
include: "rules/contacts.smk"
include: "rules/read_analysis.smk"
include: "rules/stats.smk"