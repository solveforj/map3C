
rule filter_star_pe:
    input:
        rules.remove_dna_reads_pe.output.bam
    output:
        temp("{id}_STAR_PE_filtered.bam")
    threads:
        10
    params:
        extra=config["rna"]["postprocess"]["filter_pe_params"]
    conda:
        "map3C_utils"
    shell:
        """
        samtools view {params.extra} -b {input} > {output}
        """

rule filter_star_se1:
    input:
        rules.remove_dna_reads_se1.output.bam
    output:
        temp("{id}_STAR_SE1_filtered.bam")
    threads:
        10
    params:
        extra=config["rna"]["postprocess"]["filter_se_params"]
    conda:
        "map3C_utils"
    shell:
        """
        samtools view {params.extra} -b {input} > {output}
        """

rule filter_star_se2:
    input:
        rules.remove_dna_reads_se2.output.bam
    output:
        temp("{id}_STAR_SE2_filtered.bam")
    threads:
        10
    params:
        extra=config["rna"]["postprocess"]["filter_se_params"]
    conda:
        "map3C_utils"
    shell:
        """
        samtools view {params.extra} -b {input} > {output}
        """

# Remove duplicates 

rule rmdup_star_pe:
    input:
        rules.filter_star_pe.output
    output:
        bam="{id}_STAR_PE_dedup.bam",
        stats=temp("{id}_STAR_PE_picard_log.txt")
    threads:
        10
    params:
        extra=config["rna"]["postprocess"]["picard_params"]
    conda:
        "map3C_rna"
    shell:
        """
        picard MarkDuplicates -I {input} \
            --ASSUME_SORT_ORDER "queryname" --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --ADD_PG_TAG_TO_READS false --REMOVE_DUPLICATES \
            -O {output.bam} -M {output.stats}
        """

rule rmdup_star_se1:
    input:
        rules.filter_star_se1.output
    output:
        bam="{id}_STAR_SE1_dedup.bam",
        stats=temp("{id}_STAR_SE1_picard_log.txt")
    threads:
        10
    params:
        extra=config["rna"]["postprocess"]["picard_params"]
    conda:
        "map3C_rna"
    shell:
        """
        picard MarkDuplicates -I {input} \
            --ASSUME_SORT_ORDER "queryname" --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --ADD_PG_TAG_TO_READS false --REMOVE_DUPLICATES \
            -O {output.bam} -M {output.stats}
        """

rule rmdup_star_se2:
    input:
        rules.filter_star_se2.output
    output:
        bam="{id}_STAR_SE2_dedup.bam",
        stats=temp("{id}_STAR_SE2_picard_log.txt")
    threads:
        10
    params:
        extra=config["rna"]["postprocess"]["picard_params"]
    conda:
        "map3C_rna"
    shell:
        """
        picard MarkDuplicates -I {input} \
            --ASSUME_SORT_ORDER "queryname" --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --ADD_PG_TAG_TO_READS false --REMOVE_DUPLICATES \
            -O {output.bam} -M {output.stats}
        """

def process_rna_dedup_stats(pe, se1, se2):
    prefixes = ["Picard_PE_", "Picard_SE1_", "Picard_SE2_"]
    files = [pe, se1, se2]
    out_cols = []
    out_vals = []
    for i in range(len(files)):
        file = files[i]
        prefix = prefixes[i]
        with open(file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                if line.startswith("LIBRARY"):
                    cols = [prefix + i.lower() for i in line.strip().split("\t")][1:]
                    vals = f.readline()
                    vals = vals.strip().split("\t")[1:]
                    if i > 0:
                        cols = cols[:-1]
                    out_cols += cols
                    out_vals += vals
                    break
    return out_cols, out_vals

rule process_picard_rna_dedup_stats:
    input:
        pe=rules.rmdup_star_pe.output.stats,
        se1=rules.rmdup_star_se1.output.stats,
        se2=rules.rmdup_star_se2.output.stats
    output:
        stats=temp("{id}_picard_rna_stats.txt")
    threads:
        1
    run:
        all_cols, all_vals = process_rna_dedup_stats(input["pe"], input["se1"], input["se2"])
        with open(output["stats"], "w") as f:
            f.write("\t".join(all_cols) + "\n")
            f.write("\t".join(all_vals) + "\n")        

# Gene counts

rule gene_count_pe:
    input:
        rules.rmdup_star_pe.output.bam
    output:
        counts=temp("{id}_STAR_PE_gene_counts.txt"),
        stats=temp("{id}_STAR_PE_gene_counts.txt.summary")
    threads:
        10
    params:
        extra=config["rna"]["postprocess"]["featurecounts_pe_params"]
    conda:
        "map3C_rna"
    shell:
        """
        featureCounts {params.extra} -T {threads} -t gene -o {output.counts} --donotsort {input}
        """

rule gene_count_se1:
    input:
        rules.rmdup_star_se1.output.bam
    output:
        counts=temp("{id}_STAR_SE1_gene_counts.txt"),
        stats=temp("{id}_STAR_SE1_gene_counts.txt.summary")
    threads:
        10
    params:
        extra=config["rna"]["postprocess"]["featurecounts_se_params"]
    conda:
        "map3C_rna"
    shell:
        """
        featureCounts {params.extra} -T {threads} -t gene -o {output.counts} --donotsort {input}
        """

rule gene_count_se2:
    input:
        rules.rmdup_star_se2.output.bam
    output:
        counts=temp("{id}_STAR_SE2_gene_counts.txt"),
        stats=temp("{id}_STAR_SE2_gene_counts.txt.summary")
    threads:
        10
    params:
        extra=config["rna"]["postprocess"]["featurecounts_se_params"]
    conda:
        "map3C_rna"
    shell:
        """
        featureCounts {params.extra} -T {threads} -t gene -o {output.counts} --donotsort {input}
        """

# Exon counts

rule exon_count_pe:
    input:
        rules.rmdup_star_pe.output.bam
    output:
        counts=temp("{id}_STAR_PE_exon_counts.txt"),
        stats=temp("{id}_STAR_PE_exon_counts.txt.summary")
    threads:
        10
    params:
        extra=config["rna"]["postprocess"]["featurecounts_pe_params"]
    conda:
        "map3C_rna"
    shell:
        """
        featureCounts {params.extra} -T {threads} -t exon -o {output.counts} --donotsort {input}
        """

rule exon_count_se1:
    input:
        rules.rmdup_star_se1.output.bam
    output:
        counts=temp("{id}_STAR_SE1_exon_counts.txt"),
        stats=temp("{id}_STAR_SE1_exon_counts.txt.summary")
    threads:
        10
    params:
        extra=config["rna"]["postprocess"]["featurecounts_se_params"]
    conda:
        "map3C_rna"
    shell:
        """
        featureCounts {params.extra} -T {threads} -t exon -o {output.counts} --donotsort {input}
        """

rule exon_count_se2:
    input:
        rules.rmdup_star_se2.output.bam
    output:
        counts=temp("{id}_STAR_SE2_exon_counts.txt"),
        stats=temp("{id}_STAR_SE2_exon_counts.txt.summary")
    threads:
        10
    params:
        extra=config["rna"]["postprocess"]["featurecounts_se_params"]
    conda:
        "map3C_rna"
    shell:
        """
        featureCounts {params.extra} -T {threads} -t exon -o {output.counts} --donotsort {input}
        """
    
rule aggregate_exon_counts:
    input:
        pe=rules.exon_count_pe.output.counts,
        se1=rules.exon_count_se1.output.counts,
        se2=rules.exon_count_se2.output.counts
    output:
        counts="{id}_STAR_exon_counts.txt"
    threads:
        1
    run:
        with open(input["pe"]) as pe_f, \
              open(input["se1"]) as se1_f, \
              open(input["se2"]) as se2_f, \
              open(output["counts"], "w") as counts_f:

            header=["Geneid", "Chr", "Start", "End", "Strand", "Length", "PE", "SE1", "SE2"]
            counts_f.write("\t".join(header) + "\n")

            for lines in zip(pe_f, se1_f, se2_f):
                if lines[0].startswith("#") or lines[0].startswith("Geneid"):
                    continue
                pe_l = lines[0].strip().split()
                se1_l = lines[1].strip().split()
                se2_l = lines[2].strip().split()
    
                keep = pe_l + [se1_l[-1], se2_l[-1]]
                counts_f.write("\t".join(keep) + "\n")


rule aggregate_gene_counts:
    input:
        pe=rules.gene_count_pe.output.counts,
        se1=rules.gene_count_se1.output.counts,
        se2=rules.gene_count_se2.output.counts
    output:
        counts="{id}_STAR_gene_counts.txt"
    threads:
        1
    run:
        with open(input["pe"]) as pe_f, \
              open(input["se1"]) as se1_f, \
              open(input["se2"]) as se2_f, \
              open(output["counts"], "w") as counts_f:

            header=["Geneid", "Chr", "Start", "End", "Strand", "Length", "PE", "SE1", "SE2"]
            counts_f.write("\t".join(header) + "\n")

            for lines in zip(pe_f, se1_f, se2_f):
                if lines[0].startswith("#") or lines[0].startswith("Geneid"):
                    continue
                pe_l = lines[0].strip().split()
                se1_l = lines[1].strip().split()
                se2_l = lines[2].strip().split()

                if pe_l[0] != se1_l[0] != se2_l[0]:
                    raise Exception("Gene/Exon IDs don't match!")
    
                keep = pe_l + [se1_l[-1], se2_l[-1]]
                counts_f.write("\t".join(keep) + "\n")

def process_featureCount_stats(files, prefixes):
    all_cols = []
    all_vals = []
    for i in range(len(files)):
        file = files[i]
        prefix = prefixes[i]

        with open(file) as f:
            for line in f:
                if "Status" in line:
                    continue
                line = line.strip().split()
                all_cols.append(f"{prefix}{line[0].lower()}")
                all_vals.append(f"{line[1]}")
    return all_cols, all_vals

rule aggregate_featureCount_stats:
    input:
        pe_exon=rules.exon_count_pe.output.stats,
        se1_exon=rules.exon_count_se1.output.stats,
        se2_exon=rules.exon_count_se2.output.stats,
        pe_gene=rules.gene_count_pe.output.stats,
        se1_gene=rules.gene_count_se1.output.stats,
        se2_gene=rules.gene_count_se2.output.stats,        
    output:
        stats=temp("{id}_featureCount_stats.txt")
    threads:
        1
    run:
        files = [input["pe_exon"], input["se1_exon"], input["se2_exon"],
                input["pe_gene"], input["se1_gene"], input["se2_gene"]]
        prefixes = ["featureCount_PE_exon_", "featureCount_SE1_exon_", "featureCount_SE2_exon_",
                    "featureCount_PE_gene_", "featureCount_SE1_gene_", "featureCount_SE2_gene_"]

        all_cols, all_vals = process_featureCount_stats(files, prefixes)
        with open(output["stats"], "w") as f:
            f.write("\t".join(all_cols) + "\n")
            f.write("\t".join(all_vals) + "\n")   
