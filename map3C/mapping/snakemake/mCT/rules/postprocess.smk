
if dna_done:

    dna_postprocess_protocol = config["dna"]["postprocess"]["postprocess_protocol"]

    if dna_postprocess_protocol.endswith("smk"):

        include: dna_postprocess_protocol

    elif dna_postprocess_protocol == "default":
    
        include: "postprocess_dna.smk"
    
    else:

        raise Exception("DNA alignment protocol not available for mode")

if rna_done:

    rna_postprocess_protocol = config["rna"]["postprocess"]["postprocess_protocol"]

    if rna_postprocess_protocol.endswith("smk"):

        include: rna_postprocess_protocol

    elif rna_postprocess_protocol == "default":
    
        include: "postprocess_rna.smk"
    
    else:

        raise Exception("RNA alignment protocol not available for mode")

