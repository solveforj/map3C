
rule sort_contacts:
    input:
        contacts = rules.generate_contacts.output.contacts
    output:
        contacts_sorted = (
            "{id}_contacts.sorted.pairs.gz"
            if last_contacts_step == "sort"
            else temp("{id}_contacts.sorted.pairs.gz")
        )
    params:
        extra=config["contacts"]["sort"]["sort_params"],
    conda:
        "map3C_utils"
    threads:
        1
    shell:
        """
        pairtools sort {params.extra} --nproc {threads} --output {output.contacts_sorted} {input.contacts} 
        """

rule sort_artefacts:
    input:
        artefacts = rules.generate_contacts.output.artefacts
    output:
        artefacts_sorted = (
            "{id}_artefacts.sorted.pairs.gz"
            if last_contacts_step == "sort"
            else temp("{id}_artefacts.sorted.pairs.gz")
        )
    params:
        extra=config["contacts"]["sort"]["sort_params"],
    conda:
        "map3C_utils"
    threads:
        1
    shell:
        """
        pairtools sort {params.extra} --nproc {threads} --output {output.artefacts_sorted} {input.artefacts} 
        """


