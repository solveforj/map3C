# Installation

Use conda/mamba to install the following environments:

* [map3C_preprocess_cutadapt](../../envs/preprocess/map3C_preprocess_GAGE-seq.yml)
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

Make a directory for this map3C run:

```{bash}
mkdir /path/to/map3C_run
```
You will need to index your reference genomes and generate restriction enzyme site position files. 

For reference genome, you will need to run the following commands:

```{bash}
conda activate map3C_utils
bwa index /path/to/ref.fa
```

For restriction enzyme site position files, you will need to run the following commands

```{bash}
conda activate map3C_tools
map3C restriction-sites --cut-seqs TTAA GTAC --reference /path/to/ref.fa --output /path/to/map3C_run/txt/MseI_CviQI.txt
```

Also, you will need to download the chromosome size files for your reference genome.

# Running map3C (demultiplexing)

> _A note on FASTQ file formats_
> 
>FASTQ file names should follow the format of SampleName_R1.fastq.gz. Here:
> * SampleName is the sample name and it should contain no underscores (“_”)
> * R1 is the read. In this example, R1 means Read 1. For a paired-end run, there is at least one file with R2 in the file name for Read 2. When generated, index reads are I1 or I2.

1. Update [`txt/plate_info.txt`](txt/plate_info.txt)
   * TSV where first column is sample name (no underscores allowed - Illumina convention) and second column is directory where FASTQ files for the plate are stored
2. Update [`txt/demultiplex_config_GAGE-seq.yml`](txt/demultiplex_config_GAGE-seq.yml) with:
   * The path to [`txt/plate_info.txt`](txt/plate_info.txt) goes in the _general -> fastq_info_ entry
   * The path to demultiplex directory - this will include your demultiplexed FASTQ files and no mapping results - goes in the _general -> output_directory_ entry
   * GAGE-seq BC1 barcode FASTA (the proper file is provided at [`txt/barcodes_BC1.fa`](txt/barcodes_BC1.fa)) goes in the _demultiplex_protocols -> GAGE-seq -> BC1 -> barcodes_ entry
   * GAGE-seq BC2 barcode FASTA (the proper file is provided at [`txt/barcodes_BC2.fa`](txt/barcodes_BC2.fa)) goes in the _demultiplex_protocols -> GAGE-seq -> BC2 -> barcodes_ entry
3. Run [`scripts/A01a_prepare_demultiplex.sh`](scripts/A01a_prepare_demultiplex.sh)
   * Make sure to update the path to [`txt/demultiplex_config.yml`](txt/demultiplex_config_snm3C.yml)
   * This is fast
4. Run [`scripts/A01b_run_BC2.sh`](scripts/A01b_run_BC2.sh)
   * This performs the first level of demultiplexing
   * This should be submitted with qsub (for SGE, not sure command for other systems)
   * Don’t forget to make sure that you run a job array of the correct length (number of samples) - i.e. if you have 2 samples, the qsub parameter should be -t 1-2:1
   * A01b should finish in a couple hours for each sample. Check to make sure that there is a file for each sample that has this format: `demultiplex/results/{sample}/{sample}-BC2_demultiplex_stats.txt`. This file's text should indicate that the sample finished demultiplexing in a specified amount of time.
5. Run [`scripts/A01c_run_BC1.sh`](scripts/A01c_run_BC1.sh)
   * This performs the second level of demultiplexing
   * This should be submitted with qsub (for SGE, not sure command for other systems)
   * Don’t forget to make sure that you run a job array of the correct length (number of samples x 96) - i.e. if you have 2 samples, the qsub parameter should be -t 1-192:1
   * A01c should finish in a few minutes for each job in the array. Check to make sure that there is a file for each job that has this format: `demultiplex/results/{sample}/{sample}-BC2-{BC2}/{sample}-BC2-{BC2}-BC1_demultiplex_stats.txt`. This file's text should indicate that the sample finished demultiplexing in a specified amount of time.

# Running map3C (mapping)

1. Update [`txt/mapping_info.txt`](txt/mapping_info.txt)
   * TSV where first column is cell name (underscores are allowed), second column is the whole path to R1 FASTQ, and third column is the whole path to R2 FASTQ
   * Each cell is a unique combination of BC1 and BC2 barcodes. Recommended to filter for cells that have a sufficient number of reads in their FASTQ files, as these are more likely to be "real"
2. Update [`txt/mapping_config_GAGE-seq.yml`](txt/mapping_config_GAGE-seq.yml)
   * Don’t forget to specify the correct location of mapping_info.txt in the _general -> fastq_info_ entry
   * Don’t forget to specify your mapping directory (should be different from demultiplex directory) in the _general -> output_directory_ entry
   * Go to the _align -> align_params -> bwa -> reference_path_ entry and make sure proper reference genome paths are specified
   * Go to the _contacts -> call -> call_params_ entry and make sure proper chrom sizes and cut site files are specified.
3. Run [`scripts/A02a_prepare_mapping.sh`](scripts/A02a_prepare_mapping.sh)
   * Don’t forget to specify the correct location of [`txt/mapping_config_GAGE-seq.yml`](txt/mapping_config_GAGE-seq.yml)
   * This is fast
4. Run [`scripts/A02b_run_mapping.sh`](scripts/A02b_run_mapping.sh)
   * This should be submitted with qsub (for SGE, not sure command for other systems)
   * Make sure that path to mapping_scripts.txt is correct (depends on what you named your mapping directory)
   * Don’t forget to make sure that you run a job array of the correct length, which is the number of cells