
rule trim_fastq:
    input:
        unpack(get_fastq)
    output:
        r1=temp("{id}_trimmed_R1.fastq.gz"),
        r2=temp("{id}_trimmed_R2.fastq.gz"),
        stats=temp("{id}_cutadapt_stats.txt")
    threads:
        2
    conda:
        "map3C_preprocess_cutadapt"
    params:
        extra=config["trim_methods"][trim_protocol]["cutadapt_params"]
    shell:
        """
        cutadapt --report=minimal -j {threads} {params.extra} -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {output.stats}
        """
        
rule trim_stats:
    input:
        stats=rules.trim_fastq.output.stats
    output:
        stats=temp("{id}_trim_stats.txt")
    threads:
        1
    run:
        with open(input["stats"]) as f, open(output["stats"], "w") as out:
            lines = f.readlines()
            tstats = lines[1].strip().split("\t")
            in_pairs = tstats[1]
            out_pairs = tstats[6]
            out.write("\t".join(["pre_trimmed_pairs", "trimmed_pairs"]) + "\n")
            out.write("\t".join([in_pairs, out_pairs]) + "\n")

trim_output = "separate"

def get_trimmed_r1_fastq(wildcards):
    return f"{wildcards.id}_trimmed_R1.fastq.gz"

def get_trimmed_r2_fastq(wildcards):
    return f"{wildcards.id}_trimmed_R2.fastq.gz"
    