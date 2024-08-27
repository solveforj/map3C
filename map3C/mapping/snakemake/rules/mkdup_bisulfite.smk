
rule mark_duplicates:
    input:
        get_merged_bam
    output:
        bam = temp("{id}_mkdup.bam"),
        stats = temp("{id}_dupsifter_stats.tmp.txt")
    params:
        reference_path=config["align"]["align_params"]["biscuit"]["reference_path"],
        extra=config["read_duplicates"]["dupsifter_params"]
    conda:
        "map3C_utils"
    threads:
        2
    shell:
        """
        dupsifter -s -o {output.bam} -O {output.stats} {params.extra} {params.reference_path}  {input} 
        """

rule coord_sort_mkdup:
    input:
        rules.mark_duplicates.output.bam
    output:
        temp("{id}_mkdup_sorted.bam")
    params:
        out_prefix=lambda wildcards: f"{wildcards.id}"
    conda:
        "map3C_utils"
    threads: 
        10
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        """

rule duplicates_stats:
    input:
        stats = rules.mark_duplicates.output.stats
    output:
        stats = temp("{id}_dupsifter_stats.txt")
    params:
        out_prefix=lambda wildcards: f"{wildcards.id}"
    threads:
        1
    run:
        dup_count = 0
        with open(input["stats"]) as f:
            line_count = 0
            for line in f:
                if line_count in [6, 7]:
                    dup_count += int(line.strip().split(": ")[-1])
                line_count += 1
        with open(output["stats"], "w") as f:
            f.write("duplicate_mates\n")
            f.write(str(dup_count) + "\n")

def get_merged_bam(wildcards):
    return f"{wildcards.id}_mkdup.bam"
