
rule dedup_contacts:
    input:
        contacts = rules.sort_contacts.output.contacts_sorted
    output:
        contacts_dedup = (
            "{id}_contacts.dedup.pairs.gz"
            if last_contacts_step == "dedup"
            else temp("{id}_contacts.dedup.pairs.gz")
        ),
        stats = temp("{id}_contacts_dedup_stats.txt")
    params:
        extra=config["contacts"]["dedup"]["dedup_params"]
    conda:
        "map3C_utils"
    threads:
        1
    shell:
        """
        pairtools dedup {params.extra} -p {threads} --output {output.contacts_dedup} --output-stats {output.stats} {input.contacts}
        """

rule dedup_artefacts:
    input:
        artefacts = rules.sort_artefacts.output.artefacts_sorted
    output:
        artefacts_dedup = (
            "{id}_artefacts.dedup.pairs.gz"
            if last_contacts_step in ["dedup", "lowcov"]
            else temp("{id}_artefacts.dedup.pairs.gz")
        ),
        stats = temp("{id}_artefacts_dedup_stats.txt")
    params:
        extra=config["contacts"]["dedup"]["dedup_params"]
    conda:
        "map3C_utils"
    threads:
        1
    shell:
        """
        pairtools dedup {params.extra} -p {threads} --output {output.artefacts_dedup} --output-stats {output.stats} {input.artefacts}
        """

