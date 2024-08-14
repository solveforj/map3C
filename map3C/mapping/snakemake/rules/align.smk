
align_protocol = config["align"]["align_protocol"]

if align_protocol.endswith("smk"):

    include: align_protocol

elif align_protocol == "default":

    if mode == "bsdna":

        include: "align_bisulfite.smk"

    elif mode == "dna":

        include: "align_dna.smk"

    else:

        raise Exception("Alignment protocol not available for mode")

else:

    raise Exception("Alignment protocol not appropriately specified. Contacts and/or methylation cannot be called.")

    

