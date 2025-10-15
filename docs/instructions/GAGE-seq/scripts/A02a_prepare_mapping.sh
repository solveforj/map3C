
conda activate map3C_tools

map3C prepare-mapping \
    --config /u/project/jflint/jgalasso/gageseq_dev/txt/mapping_config.yml \
    --snakemake-params "--rerun-incomplete --nolock -p -c 30"

conda deactivate