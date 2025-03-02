
if config["stats"]["stats_protocol"] != "default":
    raise Exception("Stats protocol must be specified for desired output.")

def get_all_stats(wildcards):
    stats_files = []
    if config["preprocess"]["trim_protocol"] != "none":
        stats_files.append(f"{wildcards.id}_trim_stats.txt")
    if config["dna"]["run"]:
        stats_files.append(f"{wildcards.id}_biscuit_mod_contam_stats.txt")
        stats_files.append(f"{wildcards.id}_biscuit_dupsifter_stats.txt")
        stats_files.append(f"{wildcards.id}_biscuit_alignment_stats.txt")
        stats_files.append(f"{wildcards.id}_methylation_stats.txt")
    if config["rna"]["run"]:
        stats_files.append(f"{wildcards.id}_STAR_stats.txt")
        stats_files.append(f"{wildcards.id}_STAR_contam_stats.txt")
        stats_files.append(f"{wildcards.id}_picard_rna_stats.txt")
        stats_files.append(f"{wildcards.id}_featureCount_stats.txt")

    return stats_files


rule aggregate_stats: 
    input:
        get_all_stats
    output:
        stats="{id}_qc_stats.txt"
    params:
        out_prefix = lambda wildcards: f"{wildcards.id}",
        mode = mode,
        extra = config["stats"]["stats_params"]
    conda:
        "map3C_tools"
    threads:
        1
    shell:
        'map3C aggregate-qc-stats '
        '--job {params.out_prefix} '
        '--out-prefix {params.out_prefix} '
        '--mode {params.mode} '
        '{params.extra} '