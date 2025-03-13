
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

rule duplicates_stats:
    input:
        stats = rules.mark_duplicates.output.stats
    output:
        stats = temp("{id}_dupsifter_stats.txt")
    threads:
        1
    run:
        input_count = 0
        dup_count = 0
        with open(input["stats"]) as f:
            line_count = 0
            for line in f:
                if line_count == 2:
                    input_count += int(line.strip().split(": ")[-1]) * 2
                if line_count in [3, 4]:
                    input_count += int(line.strip().split(": ")[-1])
                if line_count == 5:
                    dup_count += int(line.strip().split(": ")[-1]) * 2
                if line_count in [6, 7]:
                    dup_count += int(line.strip().split(": ")[-1])
                line_count += 1
        with open(output["stats"], "w") as f:
            f.write("dupsifter_input_mapped_mates\tdupsifter_removed_duplicate_mates\n")
            f.write(f"{input_count}\t{dup_count}\n")

def get_merged_bam(wildcards):
    return f"{wildcards.id}_mkdup.bam"
