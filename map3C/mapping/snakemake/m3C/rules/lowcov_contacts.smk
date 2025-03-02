
rule pairtools_filterbycov:
    input:
        unpack(get_pairs_data)
    output:
        lowcov = (
            "{id}_map3C.srt.dedup.lcov.pairs.gz"
            if last_contacts_step == "lowcov"
            else temp("{id}_map3C.srt.dedup.lcov.pairs.gz")
        ),
        highcov = (
            "{id}_map3C.srt.dedup.hcov.pairs.gz"
            if keep_highcov
            else temp("{id}_map3C.srt.dedup.hcov.pairs.gz")        
        ),
        stats = temp("{id}_filterbycov_stats.txt")
    params:
        extra=config["contacts"]["lowcov"]["filterbycov_params"]
    conda:
        "map3C_utils"
    threads:
        1
    shell:
        """
        contact_count=`zcat {input.contacts} | awk '!/^#/{{count++}} END{{ print count+0 }}' -`
        if [ $contact_count -ge 1 ]; then
            pairtools filterbycov \
                --output {output.lowcov} \
                --output-highcov {output.highcov} \
                --output-stats {output.stats} \
                {input.contacts}
        else
            touch {output.highcov}
            touch {output.stats}
            cp {input.contacts} {output.lowcov} 
        fi
        """

