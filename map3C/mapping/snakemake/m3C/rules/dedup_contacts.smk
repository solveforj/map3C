
rule dedup_contacts:
    input:
        contacts = rules.sort_contacts.output.contacts_sorted
    output:
        contacts_dedup = (
            "{id}_map3C.srt.dedup.pairs.gz"
            if last_contacts_step == "dedup"
            else temp("{id}_map3C.srt.dedup.pairs.gz")
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
