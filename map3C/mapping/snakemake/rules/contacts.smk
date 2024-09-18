

if config["contacts"]["call"]["call_protocol"] != "default":
    raise Exception("Contact calling protocol must be specified for desired output.")

rule generate_contacts:
    input:
        get_merged_bam
    output:
        bam = (
            "{id}_trimmed.bam"
            if bam_generated
            else []        
        ),
        contacts=(
            "{id}_all.pairs.gz"
            if last_contacts_step == "call"
            else temp("{id}_all.pairs.gz")
        ),
        stats=temp("{id}_alignment_stats.txt")
    params:
        out_prefix=lambda wildcards: f"{wildcards.id}",
        extra=config["contacts"]["call"]["call_params"],
        manual_mate_annotation=('--manual-mate-annotation ' 
                                if trim_output == "separate" and not joint_alignments 
                                else ''),
        read_type="wgs" if mode == "dna" else "bisulfite"
    conda:
        "map3C_tools"
    threads:
        1
    shell:
        'map3C call-contacts '
        '--bam {input} '
        '--out-prefix {params.out_prefix} '
        '--read-type {params.read_type} '
        '{params.manual_mate_annotation} '
        '{params.extra} '


def get_pairs_data(wildcards):
   return {"contacts" : f"{wildcards.id}_all.pairs.gz"}
    
def get_pairs_stats(wildcards):
   return {"contacts_stats": []}

def get_filterbycov_stats(wildcards):
    return {"filterbycov_stats": []}

if config["contacts"]["sort"]["sort_protocol"] == "default":

    include: "sort_contacts.smk"

    def get_pairs_data(wildcards):
       return {"contacts" : f"{wildcards.id}_all.srt.pairs.gz"}

if config["contacts"]["dedup"]["dedup_protocol"] == "default":

    if config["contacts"]["sort"]["sort_protocol"] == "none":

        raise Exception("In order to dedup contacts, they must be sorted.")

    if "--no-flip" in config["contacts"]["call"]["call_params"]:

        raise Exception("In order to dedup contacts, they must be triu-flipped.")

    include: "dedup_contacts.smk"

    def get_pairs_data(wildcards):
       return {"contacts" : f"{wildcards.id}_all.srt.dedup.pairs.gz"}
        
    def get_pairs_stats(wildcards):
       return {"contacts_stats": f"{wildcards.id}_contacts_dedup_stats.txt"}

if config["contacts"]["lowcov"]["lowcov_protocol"] == "default":

    if config["contacts"]["sort"]["sort_protocol"] == "none":

        raise Exception("In order to filterbycov contacts, they must be sorted.")

    if config["contacts"]["dedup"]["dedup_protocol"] == "none":

        raise Exception("In order to filterbycov contacts, this pipeline requires that they are deduplicated.")

    if "--no-flip" in config["contacts"]["call"]["call_params"]:

        raise Exception("In order to filterbycov contacts, they must be triu-flipped.")

    include: "lowcov_contacts.smk"
    
    def get_pairs_data(wildcards):
       return {"contacts" : f"{wildcards.id}_all.srt.dedup.lcov.pairs.gz"}
        
    def get_filterbycov_stats(wildcards):
        return {"filterbycov_stats": f"{wildcards.id}_filterbycov_stats.txt"}

if config["contacts"]["filter"]["filter_protocol"] == "default":

    include: "filter_contacts.smk"
    
    def get_pairs_data(wildcards):
       return {"contacts" : f"{wildcards.id}_{filter_suffix}.flt.pairs.gz"}

rule pairtools_stats:
    input:
        unpack(get_pairs_data),
        unpack(get_pairs_stats),
        unpack(get_filterbycov_stats)
    output:
        temp("{id}_pairs_stats.txt")
    params:
        out_prefix=lambda wildcards: f"{wildcards.id}"
    threads:
        1
    conda:
        "map3C_tools"
    shell:
        "map3C pairtools-stats "
        "--out-prefix {params.out_prefix} "
        "--contacts {input.contacts} "
        "--contacts-dedup-stats {input.contacts_stats} "
        "--filterbycov-stats {input.filterbycov_stats} "
