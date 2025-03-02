

if config["contacts"]["call"]["call_protocol"] != "default":
    raise Exception("Contact calling protocol must be specified for desired output.")


generate_contacts_bam = "{id}_map3C.bam" if bam_generated else []

rule generate_contacts:
    input:
        get_merged_bam
    output:
        bam = (
            generate_contacts_bam
            if last_bam_step == "call"
            else temp(generate_contacts_bam)        
        ),
        contacts=(
            "{id}_map3C.pairs.gz"
            if last_contacts_step == "call"
            else temp("{id}_map3C.pairs.gz")
        ),
        stats=temp("{id}_alignment_stats.txt")
    params:
        out_prefix=lambda wildcards: f"{wildcards.id}",
        extra=config["contacts"]["call"]["call_params"],
        mate_annotation=('--mate-annotation qname ' 
                                if trim_output == "separate" and not joint_alignments 
                                else '--mate-annotation flag '),
    conda:
        "map3C_tools"
    threads:
        1
    shell:
        'map3C call-contacts '
        '--bam {input} '
        '--out-prefix {params.out_prefix} '
        '{params.mate_annotation} '
        '{params.extra} '

def get_namesorted_analysis_bam(wildcards):
    return f"{wildcards.id}_map3C.bam"

def get_pairs_data(wildcards):
   return {"contacts" : f"{wildcards.id}_map3C.pairs.gz"}
    
def get_pairs_stats(wildcards):
   return {"contacts_stats": []}

def get_filterbycov_stats(wildcards):
    return {"filterbycov_stats": []}

if config["contacts"]["sort"]["sort_protocol"] == "default":

    include: "sort_contacts.smk"

    def get_pairs_data(wildcards):
       return {"contacts" : f"{wildcards.id}_map3C.srt.pairs.gz"}

if config["contacts"]["dedup"]["dedup_protocol"] == "default":

    if config["contacts"]["sort"]["sort_protocol"] == "none":

        raise Exception("In order to dedup contacts, they must be sorted.")

    if "--no-flip" in config["contacts"]["call"]["call_params"]:

        raise Exception("In order to dedup contacts, they must be triu-flipped.")

    include: "dedup_contacts.smk"

    def get_pairs_data(wildcards):
       return {"contacts" : f"{wildcards.id}_map3C.srt.dedup.pairs.gz"}
        
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
       return {"contacts" : f"{wildcards.id}_map3C.srt.dedup.lcov.pairs.gz"}
        
    def get_filterbycov_stats(wildcards):
        return {"filterbycov_stats": f"{wildcards.id}_filterbycov_stats.txt"}

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
        "map3C pair-stats "
        "--out-prefix {params.out_prefix} "
        "--input-pairs {input.contacts} "
        "--pairs-dedup-stats {input.contacts_stats} "
        "--pairs-filterbycov-stats {input.filterbycov_stats} "
