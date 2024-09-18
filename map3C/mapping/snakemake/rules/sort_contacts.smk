
rule sort_contacts:
    input:
        contacts = rules.generate_contacts.output.contacts
    output:
        contacts_sorted = (
            "{id}_all.srt.pairs.gz"
            if last_contacts_step == "sort"
            else temp("{id}_all.srt.pairs.gz")
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

