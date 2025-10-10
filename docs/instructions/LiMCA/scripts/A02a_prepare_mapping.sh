
conda activate map3C_tools

map3C prepare-mapping \
    --config /path/to/map3C_run/txt/mapping_config.yml \
    --snakemake-params "--rerun-incomplete --nolock -p -c 30"

conda deactivate