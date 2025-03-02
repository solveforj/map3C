
rule align_star:
    input:
        r1=get_trimmed_r1_fastq, 
        r2=get_trimmed_r2_fastq
    output:
        bam_pe=temp("{id}_STAR_PE_Aligned.out.bam"),
        stats_pe=temp("{id}_STAR_PE_Log.final.out"),
        bam_se1=temp("{id}_STAR_SE1_Aligned.out.bam"),
        stats_se1=temp("{id}_STAR_SE1_Log.final.out"),
        bam_se2=temp("{id}_STAR_SE2_Aligned.out.bam"),
        stats_se2=temp("{id}_STAR_SE2_Log.final.out"),
    threads: 
        10
    params:
        reference_path=config["rna"]["align"]["star"]["reference_path"],
        extra=config["rna"]["align"]["star"]["star_params"],
        load_prefix=lambda wildcards: f"{wildcards.id}_STAR_LOAD_",
        pe_prefix=lambda wildcards: f"{wildcards.id}_STAR_PE_",
        se1_prefix=lambda wildcards: f"{wildcards.id}_STAR_SE1_",
        se2_prefix=lambda wildcards: f"{wildcards.id}_STAR_SE2_",
    conda:
        "map3C_rna"
    retries: 3
    shell:
        """
        STAR --outFileNamePrefix {params.load_prefix} --genomeDir {params.reference_path} --genomeLoad LoadAndExit
        
        STAR {params.extra} --genomeDir {params.reference_path} \
            --runThreadN {threads} \
            --genomeLoad LoadAndKeep \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.pe_prefix} \
            --readFilesIn {input.r1} {input.r2} \
            --outReadsUnmapped Fastx
        
        STAR {params.extra} --genomeDir {params.reference_path} \
            --runThreadN {threads} \
            --genomeLoad LoadAndKeep \
            --outFileNamePrefix {params.se1_prefix} \
            --readFilesIn {params.pe_prefix}Unmapped.out.mate1
        
        STAR {params.extra} --genomeDir {params.reference_path} \
            --runThreadN {threads} \
            --genomeLoad LoadAndKeep \
            --outFileNamePrefix {params.se2_prefix} \
            --readFilesIn {params.pe_prefix}Unmapped.out.mate2
                
        STAR --outFileNamePrefix {params.load_prefix} --genomeDir {params.reference_path} --genomeLoad Remove

        rm *STAR*Unmapped.out.mate* 
        rm -rf *_STARtmp/ 
        rm *STAR*Aligned.out.sam *STAR*Log.out *STAR*Log.progress.out *STAR*SJ.out.tab
        rm *LOAD*
        """

def process_star_stats(stats, prefixes):
    out_cols = []
    out_vals = []
    for s in range(len(stats)):
        stat = stats[s]
        prefix = prefixes[s]
        with open(stat) as f:
            count = 0
            for line in f:
                line = line.strip()
                if "|" not in line:
                    continue
                if "%" in line:
                    continue
                count += 1
                if count < 5:
                    continue
                line = [i.strip() for i in line.split("|")]
                line[-1] = str(float(line[-1]))
                line[0] = line[0].lower().replace(" ", "_")
                line[0] = line[0].replace(":", "")
                line[0] = f"{prefix}_{line[0]}"
                out_cols.append(line[0])
                out_vals.append(line[1])
    return out_cols, out_vals

rule star_stats:
    input:
        pe=rules.align_star.output.stats_pe,
        se1=rules.align_star.output.stats_se1,
        se2=rules.align_star.output.stats_se2
    output:
        stats=temp("{id}_STAR_stats.txt")
    threads:
        1
    run:
        stats = [input["pe"], input["se1"], input["se2"]]
        prefixes = ["STAR_PE", "STAR_SE1", "STAR_SE2"]
        all_cols, all_vals = process_star_stats(stats, prefixes)
        with open(output["stats"], "w") as f:
            f.write("\t".join(all_cols) + "\n")
            f.write("\t".join(all_vals) + "\n")

rule remove_dna_reads_pe:
    input:
        rules.align_star.output.bam_pe
    output:
        bam=temp("{id}_STAR_PE_contam_filtered.bam"),
        stats=temp("{id}_STAR_PE_contam_stats.txt"),
    conda:
        "map3C_tools"
    params:
        out_prefix=lambda wildcards: f"{wildcards.id}_STAR_PE",
        extra=config["rna"]["align"]["contamination_filter_pe_params"]
    shell:
        'map3C contamination-filter '
        '--bam {input} '
        '--out-prefix {params.out_prefix} '
        '--mate-annotation flag '
        '{params.extra} '

rule remove_dna_reads_se1:
    input:
        rules.align_star.output.bam_se1
    output:
        bam=temp("{id}_STAR_SE1_contam_filtered.bam"),
        stats=temp("{id}_STAR_SE1_contam_stats.txt"),
    conda:
        "map3C_tools"
    params:
        out_prefix=lambda wildcards: f"{wildcards.id}_STAR_SE1",
        extra=config["rna"]["align"]["contamination_filter_se_params"]
    shell:
        'map3C contamination-filter '
        '--bam {input} '
        '--out-prefix {params.out_prefix} '
        '--mate-annotation flag '
        '{params.extra} '

rule remove_dna_reads_se2:
    input:
        rules.align_star.output.bam_se2
    output:
        bam=temp("{id}_STAR_SE2_contam_filtered.bam"),
        stats=temp("{id}_STAR_SE2_contam_stats.txt"),
    conda:
        "map3C_tools"
    params:
        out_prefix=lambda wildcards: f"{wildcards.id}_STAR_SE2",
        extra=config["rna"]["align"]["contamination_filter_se_params"]
    shell:
        'map3C contamination-filter '
        '--bam {input} '
        '--out-prefix {params.out_prefix} '
        '--mate-annotation flag '
        '{params.extra} '

def process_rna_contam_stats(pe, se1, se2):
    out_cols = []
    out_vals = []
    with open(pe) as f:
        cols = f.readline()
        cols = ["STAR_PE_" + i for i in cols.strip().split()]
        vals = f.readline()
        vals = vals.strip().split()
    out_cols += cols
    out_vals += vals
    with open(se1) as f:
        cols = f.readline()
        cols = ["STAR_SE1_" + i for i in cols.strip().split()][0:2]
        vals = f.readline()
        vals = vals.strip().split()[0:2]
    out_cols += cols
    out_vals += vals
    with open(se2) as f:
        cols = f.readline()
        cols = ["STAR_SE2_" + i.replace("R1", "R2") for i in cols.strip().split()][0:2]
        vals = f.readline()
        vals = vals.strip().split()[0:2]
    out_cols += cols
    out_vals += vals
    return out_cols, out_vals

rule rna_contam_stats:
    input:
        pe=rules.remove_dna_reads_pe.output.stats,
        se1=rules.remove_dna_reads_se1.output.stats,
        se2=rules.remove_dna_reads_se2.output.stats
    output:
        stats=temp("{id}_STAR_contam_stats.txt")
    threads:
        1
    run:
        all_cols, all_vals = process_rna_contam_stats(input["pe"], input["se1"], input["se2"])
        with open(output["stats"], "w") as f:
            f.write("\t".join(all_cols) + "\n")
            f.write("\t".join(all_vals) + "\n")

