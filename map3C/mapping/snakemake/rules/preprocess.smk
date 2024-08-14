
def get_fastq(wildcards):
    fq = run_info.loc[wildcards.id]
    out = {"r1" : fq["r1"], "r2" : fq["r2"]}
    return out

trim_protocol = config["preprocess"]["trim_protocol"]

if trim_protocol.endswith("smk"):

    include: trim_protocol

elif trim_protocol == "snm3Cseq":

    include: "trim/trim_snm3Cseq.smk"

elif trim_protocol == "meta":

    include: "trim/trim_meta.smk"

elif trim_protocol == "hires":

    include: "trim/trim_hires.smk"

elif trim_protocol == "none":

    trim_output = "separate"

    def get_trimmed_r1_fastq(wildcards):
        return run_info.loc[wildcards.id]["r1"]

    def get_trimmed_r2_fastq(wildcards):
        return run_info.loc[wildcards.id]["r2"]

else:

    raise Exception("Trim protocol not appropriately specified")

