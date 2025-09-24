
conda activate map3C_tools

map3C prepare-demultiplex \
    --config /path/to/map3C_run/txt/demultiplex_config_snm3C.yml \
    --snakemake-params "--rerun-incomplete --nolock -p -c 30"

conda deactivate