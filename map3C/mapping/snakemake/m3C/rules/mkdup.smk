
duplicate_protocol = config["read_duplicates"]["duplicate_protocol"]

if duplicate_protocol.endswith("smk"):

    include: duplicate_protocol

elif duplicate_protocol == "default":

    if mode == "bsdna":

        include: "mkdup_bisulfite.smk"


