
if dna_done:

    dna_align_protocol = config["dna"]["align"]["align_protocol"]

    if dna_align_protocol.endswith("smk"):

        include: dna_align_protocol

    elif dna_align_protocol == "default":
    
        include: "align_biscuit.smk"
    
    else:

        raise Exception("DNA alignment protocol not available for mode")

if rna_done:

    rna_align_protocol = config["rna"]["align"]["align_protocol"]

    if rna_align_protocol.endswith("smk"):

        include: rna_align_protocol

    elif rna_align_protocol == "default":
    
        include: "align_star.smk"
    
    else:

        raise Exception("RNA alignment protocol not available for mode")

