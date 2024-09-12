

if config["contacts"]["call"]["call_protocol"] != "default":
    raise Exception("Contact calling protocol must be specified for desired output.")

rule generate_contacts:
    input:
        get_merged_bam
    output:
        bam = (
            "{id}_trimmed.bam"
            if config["contacts"]["call"]["keep_trimmed_bam"]
            else temp("{id}_trimmed.bam")        
        ),
        contacts=(
            "{id}_contacts.pairs.gz"
            if last_contacts_step == "call"
            else temp("{id}_contacts.pairs.gz")
        ),
        artefacts=(
            "{id}_artefacts.pairs.gz"
            if last_contacts_step == "call"
            else temp("{id}_artefacts.pairs.gz")
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
   return {"contacts" : f"{wildcards.id}_contacts.pairs.gz",
           "artefacts" : f"{wildcards.id}_artefacts.pairs.gz"}
    
def get_pairs_stats(wildcards):
   return {"contacts_stats": [],
           "artefacts_stats": []}

def get_filterbycov_stats(wildcards):
    return {"filterbycov_stats": []}

def get_trimmed_bam(wildcards):
    return f"{wildcards.id}_trimmed.bam"

if config["contacts"]["sort"]["sort_protocol"] == "default":

    include: "sort_contacts.smk"

    def get_pairs_data(wildcards):
       return {"contacts" : f"{wildcards.id}_contacts.sorted.pairs.gz",
               "artefacts" : f"{wildcards.id}_artefacts.sorted.pairs.gz"
              }

if config["contacts"]["dedup"]["dedup_protocol"] == "default":

    if config["contacts"]["sort"]["sort_protocol"] == "none":

        raise Exception("In order to dedup contacts/artefacts, they must be sorted.")

    if "--no-flip" in config["contacts"]["call"]["call_params"]:

        raise Exception("In order to dedup contacts/artefacts, they must be triu-flipped.")

    include: "dedup_contacts.smk"

    def get_pairs_data(wildcards):
       return {"contacts" : f"{wildcards.id}_contacts.dedup.pairs.gz",
               "artefacts" : f"{wildcards.id}_artefacts.dedup.pairs.gz"
              }
        
    def get_pairs_stats(wildcards):
       return {"contacts_stats": f"{wildcards.id}_contacts_dedup_stats.txt",
               "artefacts_stats": f"{wildcards.id}_artefacts_dedup_stats.txt",
              }

if config["contacts"]["lowcov"]["lowcov_protocol"] == "default":

    if config["contacts"]["sort"]["sort_protocol"] == "none":

        raise Exception("In order to filterbycov contacts/artefacts, they must be sorted.")

    if config["contacts"]["dedup"]["dedup_protocol"] == "none":

        raise Exception("In order to filterbycov contacts/artefacts, this pipeline requires that they are deduplicated.")

    if "--no-flip" in config["contacts"]["call"]["call_params"]:

        raise Exception("In order to filterbycov contacts/artefacts, they must be triu-flipped.")

    include: "lowcov_contacts.smk"
    
    def get_pairs_data(wildcards):
       return {"contacts" : f"{wildcards.id}_contacts.dedup.lowcov.pairs.gz",
               "artefacts" : f"{wildcards.id}_artefacts.dedup.pairs.gz"
              }
        
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
        "map3C pairtools-stats "
        "--out-prefix {params.out_prefix} "
        "--contacts {input.contacts} "
        "--artefacts {input.artefacts} "
        "--contacts-dedup-stats {input.contacts_stats} "
        "--artefacts-dedup-stats {input.artefacts_stats} "
        "--filterbycov-stats {input.filterbycov_stats} "

rule coord_sort_trimmed:
    input:
        rules.generate_contacts.output.bam
    output:
        temp("{id}_trimmed_sorted.bam")
    threads: 
        10
    conda:
        "map3C_utils"
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        """

def get_coordsorted_analysis_bam(wildcards):
    return f"{wildcards.id}_trimmed_sorted.bam"
