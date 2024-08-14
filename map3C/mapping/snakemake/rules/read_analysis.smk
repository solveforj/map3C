
mask_protocol = config["read_analysis"]["mask"]["mask_protocol"]

if mask_protocol == "default":

    include: "mask.smk"

allc_protocol = config["read_analysis"]["allc"]["allc_protocol"]

if allc_protocol == "default":

    if mode == "bsdna":

        include: "allc.smk"




