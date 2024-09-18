
rule filter_pairs:
    input:
        unpack(get_pairs_data)
    output:
        pairs = pairs,
        split_reads = split_reads
    params:
        extra=config["contacts"]["filter"]["filter_params"],
        out_prefix=lambda wildcards: f"{wildcards.id}_{filter_suffix}"
    conda:
        "map3C_tools"
    threads:
        1
    shell:
        'map3C filter-pairs '
        '--input-pairs {input.contacts} '
        '--out-prefix {params.out_prefix} '
        '{params.extra} '

