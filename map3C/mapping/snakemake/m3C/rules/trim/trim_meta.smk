
rule trim:
    input:
        unpack(get_fastq)
    output:
        temp("{id}_trimmed.fastq.gz")
    threads:
        10
    conda:
        "map3C_preprocess_meta"
    params:
        pre_meta=config["trim_methods"]["meta"]["pre-meta"],
        extra=config["trim_methods"]["meta"]["pre-meta_params"]
    shell:
        """
        seqtk mergepe {input.r1} {input.r2} | {params.pre_meta} -t {threads} {params.extra} - | gzip > {output}
        """

rule seqtk_stats:
    input:
        rules.trim.output
    output:
        stats=temp("{id}_seqtk_stats.txt")
    threads:
        2
    conda:
        "map3C_preprocess_meta"
    shell:
        """
        seqtk size {input} > {output.stats}
        """

rule trim_stats:
    input:
        seqtk_stats=rules.seqtk_stats.output
    output:
        stats = temp("{id}_trim_stats.txt")
    threads:
        1
    run:
        with open(output["stats"], "w") as f, \
            open(input[0]) as seqtk:

            input_data = seqtk.readlines()[0].strip().split()
            out_reads = input_data[0]

            f.write("\t".join(["post_trim_reads"]) + "\n")
            f.write("\t".join([out_reads]) + "\n")
        

trim_output = "interleaved"

def get_trimmed_interleaved_fastq(wildcards):
    return f"{wildcards.id}_trimmed.fastq.gz"
