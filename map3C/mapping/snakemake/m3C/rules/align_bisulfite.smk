
contamination_protocol = config["contamination"]["contamination_protocol"]

if trim_output == "interleaved":

    rule align:
        input:
            get_trimmed_interleaved_fastq
        output:
            temp("{id}_alignments.bam")
        threads: 
            config["align"]["align_params"]["biscuit"]["threads"]
        params:
            reference_path=config["align"]["align_params"]["biscuit"]["reference_path"],
            extra=config["align"]["align_params"]["biscuit"]["joint_params"]
        conda:
            "map3C_utils"
        retries: 3
        shell:
            """
            biscuit align -SPM -@ {threads} -p {params.extra}  {params.reference_path} {input} \
                | samtools sort -o {output} -O BAM 
            """

elif trim_output == "separate" and not joint_alignments:

    rule align_r1:
        input:
            rules.reformat_r1.output
        output:
            temp("{id}_R1_alignments.bam")
        threads: 
            config["align"]["align_params"]["biscuit"]["threads"]
        params:
            reference_path=config["align"]["align_params"]["biscuit"]["reference_path"],
            extra=config["align"]["align_params"]["biscuit"]["separate_R1_params"]
        conda:
            "map3C_utils"
        retries: 3
        shell:
            """
            biscuit align -SPM -@ {threads} {params.extra}  {params.reference_path} {input} \
                | samtools sort -o {output} -O BAM 
            """
    
    rule align_r2:
        input:
            rules.reformat_r2.output
        output:
            temp("{id}_R2_alignments.bam")
        threads: 
            config["align"]["align_params"]["biscuit"]["threads"]
        params:
            reference_path=config["align"]["align_params"]["biscuit"]["reference_path"],
            extra=config["align"]["align_params"]["biscuit"]["separate_R2_params"]
        conda:
            "map3C_utils"
        retries: 3
        shell:
            """
            biscuit align -SPM -@ {threads} {params.extra} {params.reference_path} {input} \
                | samtools sort -o {output} -O BAM 
            """
    
    
    if contamination_protocol == "default":
    
        rule bsconv_r1:
            input:
                rules.align_r1.output
            output:
                temp("{id}_R1_bsconv_alignments.bam")
            threads:
                1
            params:
                reference_path=config["align"]["align_params"]["biscuit"]["reference_path"]
            conda:
                "map3C_utils"
            shell:
                """
                biscuit bsconv {params.reference_path} {input} {output} 
                """
    
        rule bsconv_r2:
            input:
                rules.align_r2.output
            output:
                temp("{id}_R2_bsconv_alignments.bam")
            threads:
                1
            params:
                reference_path=config["align"]["align_params"]["biscuit"]["reference_path"]
            conda:
                "map3C_utils"
            shell:
                """
                biscuit bsconv {params.reference_path} {input} {output} 
                """
        
        def get_bams(wildcards):
            return expand("{id}_R{mate}_bsconv_alignments.bam", id=wildcards.id, mate=[1, 2])
    
    else:
    
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
            config["align"]["align_params"]["biscuit"]["threads"]
        params:
            reference_path=config["align"]["align_params"]["biscuit"]["reference_path"],
            extra=config["align"]["align_params"]["biscuit"]["joint_params"]
        conda:
            "map3C_utils"
        retries: 3
        shell:
            """
            biscuit align -SPM -@ {threads} {params.extra} {params.reference_path} {input} \
                | samtools sort -o {output} -O BAM 
            """

if (trim_output == "interleaved") or (trim_output == "separate" and joint_alignments):

    if contamination_protocol == "default":

        rule bsconv:
            input:
                rules.align.output
            output:
                temp("{id}_bsconv_alignments.bam")
            threads:
                1
            params:
                reference_path=config["align"]["align_params"]["biscuit"]["reference_path"]
            conda:
                "map3C_utils"
            shell:
                """
                biscuit bsconv {params.reference_path} {input} {output} 
                """
    
        def get_bams(wildcards):
            return f"{wildcards.id}_bsconv_alignments.bam"

    else:
    
        def get_bams(wildcards):
            return f"{wildcards.id}_alignments.bam"
