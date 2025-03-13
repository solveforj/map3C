
rule rmdup_dna:
    input:
        rules.remove_rna_reads.output.bam
    output:
        bam=temp("{id}_biscuit_rmdup.bam"),
        stats=temp("{id}_biscuit_dup_stats.txt")
    threads: 
        10
    params:
        reference_path=config["dna"]["align"]["biscuit"]["reference_path"],
        extra=config["dna"]["postprocess"]["dupsifter_params"]
    conda:
        "map3C_utils"
    shell:
        """
        dupsifter {params.extra} -o {output.bam} -O {output.stats} {params.reference_path} {input}
        """

rule dupsifter_stats:
    input:
        stats = rules.rmdup_dna.output.stats
    output:
        stats = temp("{id}_biscuit_dupsifter_stats.txt")
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

rule generate_contacts:
    input:
        rules.rmdup_dna.output.bam
    output:
        bam = temp("{id}_biscuit_map3C.bam"),
        contacts = "{id}_biscuit_map3C.pairs.gz",
        stats=temp("{id}_biscuit_alignment_stats.txt")
    params:
        out_prefix=lambda wildcards: f"{wildcards.id}_biscuit",
        extra=config["dna"]["postprocess"]["call_contacts_params"],
    conda:
        "map3C_tools"
    threads:
        1
    shell:
        """
        map3C call-contacts \
            --bam {input} \
            --out-prefix {params.out_prefix} \
            --mate-annotation flag \
            {params.extra}
        """


rule mask:
    input:
        rules.generate_contacts.output.bam
    output:
        temp("{id}_biscuit_map3C_masked.bam")
    params:
        out_prefix=lambda wildcards: f"{wildcards.id}_biscuit",
        extra=config["dna"]["postprocess"]["mask_params"]
    conda:
        "map3C_tools"
    threads:
        1
    shell:
        """
        map3C mask-overlaps \
            --bam {input} \
            --out-prefix {params.out_prefix} \
            --mate-annotation flag \
            {params.extra}
        """


rule coord_sort_dna:
    input:
        rules.mask.output
    output:
        "{id}_biscuit_csort.bam"
    threads:
        10
    conda:
        "map3C_utils"
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        """

rule bam_to_allc:
    input:
        rules.coord_sort_dna.output
    output:
        allc = "{id}.allc.tsv.gz",
        tbi = "{id}.allc.tsv.gz.tbi",
        stats = "{id}.allc.tsv.gz.count.csv",
        methylation_stats = temp("{id}_methylation_stats.txt")
    params:
        out_prefix=lambda wildcards: f"{wildcards.id}",
        reference_path=config["dna"]["align"]["biscuit"]["reference_path"],
        extra=config["dna"]["postprocess"]["allc_params"]
    conda:
        "map3C_tools"
    threads:
        1
    shell:
        'map3C bam-to-allc '
        '--bam-path {input} '
        '--reference-fasta {params.reference_path} '
        '--out-prefix {params.out_prefix} ' 
        '--save-count-df ' 
        '{params.extra} '
