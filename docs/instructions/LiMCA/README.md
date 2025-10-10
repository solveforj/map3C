# Installation

Use conda/mamba to install the following environments:

* [map3C_preprocess_meta](../../envs/preprocess/map3C_preprocess_meta.yml)
* [map3C_snakemake](../../envs/map3C_snakemake.yml)
* [map3C_tools](../../envs/map3C_tools.yml)
* [map3C_utils](../../envs/map3C_utils.yml)

Run the following command to install map3C:

```{bash}
conda activate map3C_tools
# Feel free to specify map3C version
pip install map3C
```

# Preparation

You will need to index your reference genomes and generate restriction enzyme site position files. 

For reference genome, you will need to run the following commands:

```{bash}
conda activate map3C_utils
bwa index /path/to/ref.fa
```

For restriction enzyme site position files, you will need to run the following commands

```{bash}
conda activate map3C_tools
map3C restriction-sites --cut-seqs CATG --reference /path/to/ref.fa --output /path/to/map3C_run/txt/NlaIII.txt
```

Also, you will need to download the chromosome size files for your reference genome.

You will need to install the pre-meta software by downloading and making the executable:

```{bash}
git clone https://github.com/lh3/pre-pe.git
cd pre-pe
make
```

# Running map3C (mapping)

> _A note on FASTQ file formats_
> 
>FASTQ file names should follow the format of SampleName_R1.fastq.gz. Here:
> * SampleName is the sample name and it should contain no underscores (“_”)
> * R1 is the read. In this example, R1 means Read 1. For a paired-end run, there is at least one file with R2 in the file name for Read 2. When generated, index reads are I1 or I2.

1. Update [`txt/mapping_info.txt`](txt/mapping_info.txt)
   * TSV where first column is well name (underscores are allowed), second column is the whole path to R1 FASTQ, and third column is the whole path to R2 FASTQ
   * Should have # file lines = # of cells
2. Update [`txt/mapping_config_LiMCA.yml`](txt/mapping_config_LiMCA.yml)
   * Don’t forget to specify the correct location of mapping_info.txt in the fastq_info entry
   * Don’t forget to specify your mapping directory (should be different from demultiplex directory)
   * Go to the align section and make sure proper reference genome paths are specified
   * For LiMCA, we use BWA MEM, so this is the only one you need to change
   * Go to the contacts section and make sure proper chrom sizes and cut site files are specified.
   * Note that the _contacts -> call -> call_params_ section cut site parameters need to be in a specific order. The “restriction-sites” parameter can be specified once for each enzyme used (with a cut site location file), but the “restriction-enzymes” parameter must have the same order and must be specified once for each enzyme used.
   * Add the full path to the pre-meta executable to the _trim_methods -> meta -> pre-meta_ section
3. Run [`scripts/A02a_prepare_mapping.sh`](`scripts/A02a_prepare_mapping.sh`)
   * Don’t forget to specify the correct location of [`txt/mapping_config_HiRES.yml`](txt/mapping_config_HiRES.yml)
   * This is fast
4. Run [`scripts/A02b_run_mapping.sh`](scripts/A02b_run_mapping.sh)
   * This should be submitted with qsub (for SGE, not sure command for other systems)
   * Make sure that path to mapping_scripts.txt is correct (depends on what you named your mapping directory)
   * Don’t forget to make sure that you run a job array of the correct length - (384)x(# of plates)
