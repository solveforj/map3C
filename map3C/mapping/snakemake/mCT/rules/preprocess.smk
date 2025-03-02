
def get_fastq(wildcards):
    fq = run_info.loc[wildcards.id]
    out = {"r1" : fq["r1"], "r2" : fq["r2"]}
    return out

trim_protocol = config["preprocess"]["trim_protocol"]

if trim_protocol.endswith("smk"):

    include: trim_protocol

elif trim_protocol == "snmCTseq":

    include: "trim/trim_cutadapt_pe.smk"

elif trim_protocol == "none":

    def get_trimmed_r1_fastq(wildcards):
        return run_info.loc[wildcards.id]["r1"]

    def get_trimmed_r2_fastq(wildcards):
        return run_info.loc[wildcards.id]["r2"]

else:

    raise Exception("Trim protocol not appropriately specified")

