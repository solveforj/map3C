
joint_alignments = config["align"]["allow_joint_alignments"]   

if trim_output == "separate":

    if not joint_alignments:

        rule reformat_r1:
            input:
                get_trimmed_r1_fastq
            output:
                temp("{id}_R1_trimmed_reformat.fastq.gz"),
            threads:
                1
            run:
                if input[0].endswith(".gz"):
                    shell(
                        "gzip -cd {input} | awk -F' ' '{{l=l+1; if ((l-1)%4==0) {{$2=1;split($1,a,\"/\");print a[1] \"_\" $2;}}else{{print $0}}}}' | "
                        "gzip -c > {output}"
                    )
                else:
                    shell(
                        "cat {input} | awk -F' ' '{{l=l+1; if ((l-1)%4==0) {{$2=1;split($1,a,\"/\");print a[1] \"_\" $2;}}else{{print $0}}}}' | "
                        "gzip -c > {output}"
                    )
        
        rule reformat_r2:
            input:
                get_trimmed_r2_fastq
            output:
                temp("{id}_R2_trimmed_reformat.fastq.gz"),
            threads:
                1
            run:
                if input[0].endswith(".gz"):
                    shell(
                        "gzip -cd {input} | awk -F' ' '{{l=l+1; if ((l-1)%4==0) {{$2=2;split($1,a,\"/\");print a[1] \"_\" $2;}}else{{print $0}}}}' | "
                        "gzip -c > {output}" 
                    )
                else:
                    shell(
                        "cat {input} | awk -F' ' '{{l=l+1; if ((l-1)%4==0) {{$2=2;split($1,a,\"/\");print a[1] \"_\" $2;}}else{{print $0}}}}' | "
                        "gzip -c > {output}" 
                    )
