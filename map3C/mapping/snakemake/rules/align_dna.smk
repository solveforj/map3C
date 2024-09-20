
if trim_output == "interleaved":

    rule align:
        input:
            get_trimmed_interleaved_fastq
        output:
            temp("{id}_alignments.bam")
        threads: 
            config["align"]["align_params"]["bwa"]["threads"]
        params:
            reference_path=config["align"]["align_params"]["bwa"]["reference_path"],
            extra=config["align"]["align_params"]["bwa"]["joint_params"]
        conda:
            "map3C_utils"
        retries: 3
        shell:
            """
            bwa mem -5SPM -t {threads} -p {params.extra} {params.reference_path} {input} | samtools sort -o {output} -O BAM
            """

    def get_bams(wildcards):
    
        return f"{wildcards.id}_alignments.bam"

elif trim_output == "separate" and not joint_alignments:

    rule align_r1:
        input:
            rules.reformat_r1.output
        output:
            temp("{id}_R1_alignments.bam")
        threads: 
            config["align"]["align_params"]["bwa"]["threads"]
        params:
            reference_path=config["align"]["align_params"]["bwa"]["reference_path"],
            extra=config["align"]["align_params"]["bwa"]["separate_R1_params"]
        conda:
            "map3C_utils"
        retries: 3
        shell:
            """
            bwa mem -5SPM -t {threads} {params.extra} {params.reference_path} {input} | samtools sort -o {output} -O BAM
            """
    
    rule align_r2:
        input:
            rules.reformat_r2.output
        output:
            temp("{id}_R2_alignments.bam")
        threads: 
            config["align"]["align_params"]["bwa"]["threads"]
        params:
            reference_path=config["align"]["align_params"]["bwa"]["reference_path"],
            extra=config["align"]["align_params"]["bwa"]["separate_R2_params"]
        conda:
            "map3C_utils"
        retries: 3
        shell:
            """
            bwa mem -5SPM -t {threads} {params.extra} {params.reference_path} {input} | samtools sort -o {output} -O BAM
            """

    def get_bams(wildcards):
    
        return expand("{id}_R{mate}_alignments.bam", id=wildcards.id, mate=[1, 2])

elif trim_output == "separate" and joint_alignments:

    rule align:
        input:
            get_trimmed_r1_fastq, 
            get_trimmed_r2_fastq
        output:
            temp("{id}_alignments.bam")
        threads: 
            config["align"]["align_params"]["bwa"]["threads"]
        params:
            reference_path=config["align"]["align_params"]["bwa"]["reference_path"],
            extra=config["align"]["align_params"]["bwa"]["joint_params"]
        conda:
            "map3C_utils"
        retries: 3
        shell:
            """
            bwa mem -5SPM -t {threads} {params.extra} {params.reference_path} {input} | samtools sort -o {output} -O BAM
            """

    def get_bams(wildcards):
    
        return f"{wildcards.id}_alignments.bam"
