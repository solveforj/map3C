
rule interleave:
    input:
        unpack(get_fastq)
    output:
        temp("{id}_interleaved.fastq.gz")
    threads:
        1
    conda:
        "map3C_trim_meta"
    shell:
        """
        seqtk mergepe {input.r1} {input.r2} | gzip > {output}
        """

rule trim:
    input:
        rules.interleave.output
    output:
        temp("{id}_trimmed.fastq.gz")
    threads: 
        2
    params:
        pre_meta=config["trim_methods"]["meta"]["pre-meta"],
        extra=config["trim_methods"]["meta"]["pre-meta_params"]
    conda:
        "map3C_trim_meta"
    shell:
        """
        zcat {input} | {params.pre_meta} -t {threads} {params.extra} - | gzip > {output}
        """

rule seqtk_stats:
    input:
        rules.trim.output
    output:
        stats=temp("{id}_seqtk_stats.txt")
    threads:
        2
    conda:
        "map3C_trim_meta"
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
