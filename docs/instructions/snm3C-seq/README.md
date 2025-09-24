# Installation

Use conda/mamba to install the following environments:

* [map3C_preprocess_cutadapt](../../envs/preprocess/map3C_preprocess_cutadapt.yml)
* [map3C_snakemake](../../envs/map3C_snakemake.yml)
* [map3C_tools](../../envs/map3C_tools.yml)
* [map3C_utils](../../envs/map3C_utils.yml)

# Preparation

Make a directory for this map3C run:

```{bash}
mkdir /path/to/map3C_run
```
You will need to index your reference genomes and generate restriction enzyme site position files. 

For reference genome, you will need to run the following commands:

```{bash}
conda activate map3C_utils
biscuit index /path/to/ref.fa
```

For restriction enzyme site position files, you will need to run the following commands

```{bash}
conda activate map3C_tools
map3C restriction-sites --cut-seqs GATC --reference /path/to/ref.fa --output /path/to/map3C_run/txt/MboI.txt
map3C restriction-sites --cut-seqs CATG --reference /path/to/ref.fa --output /path/to/map3C_run/txt/NlaIII.txt
```

Also, you will need to download the chromosome size files for your reference genome.


# Running map3C (demultiplexing)

> _A note on FASTQ file formats_
> 
>The FASTQ files to be demultiplexed each contain data from a 384-well plate. Each plate has R1 and R2 files, potentially across multiple lanes. The file names should follow the Illumina convention of SampleName_S1_L001_R1_001.fastq.gz. Here:
> * SampleName is the plate name and it should contain no underscores (“_”)
> * S1 is the sample number based on the order that samples are listed in the sample sheet starting with 1. In this example, S1 indicates that this sample is the first sample listed in the sample sheet.
> * L001 is the lane number.
> * R1 is the read. In this example, R1 means Read 1. For a paired-end run, there is at least one file with R2 in the file name for Read 2. When generated, index reads are I1 or I2.
> * 001 is the last segment, which is always 001.
> 
> map3C will still work if the sample number or trailing “001” are not present, but it is imperative that the SampleName, Lane, and Read are specified.


1. Update [`txt/plate_info.txt`](txt/plate_info.txt)
   * TSV where first column is plate sample name (no underscores allowed - Illumina convention) and second column is directory where FASTQ files for the plate are stored
2. Update [`txt/demultiplex_config_snm3C.yml`](txt/demultiplex_config_snm3C.yml) with:
   * The path to [`txt/plate_info.txt`](txt/plate_info.txt) goes in the fastq_info entry
   * The path to demultiplex directory - this will include your demultiplexed FASTQ files and no mapping results - goes in the output_directory entry
   * snm3Cseq barcodes FASTA (the proper file is provided at [`txt/random_index_v2.multiplex.fa`](txt/random_index_v2.multiplex.fa))
3. Run [`scripts/A01a_prepare_demultiplex.sh`](scripts/A01a_prepare_demultiplex.sh)
   * Make sure to update the path to [`txt/demultiplex_config.yml`](txt/demultiplex_config_snm3C.yml)
   * This is fast
4. Run [`scripts/A01b_run_demultiplex.sh`](scripts/A01b_run_demultiplex.sh)
   * This should be submitted with qsub (for SGE, not sure command for other systems)
   * Don’t forget to make sure that you run a job array of the correct length (number of plates) - i.e. if you have 24 plates, the qsub parameter should be -t 1-24:1
   * A01b should finish in a couple hours for each plate. Check to make sure that there is a file for each plate that has this format: `demultiplex/results/{plate}/{plate}_demultiplex_stats.txt`

# Running map3C (mapping)

1. Update [`txt/mapping_info.txt`](txt/mapping_info.txt)
   * TSV where first column is well name (underscores are allowed), second column is the whole path to R1 FASTQ, and third column is the whole path to R2 FASTQ
   * Should have (384)x(# of plates) lines
2. Update [`txt/mapping_config_snm3C.yml`](txt/mapping_config_snm3C.yml)
   * Don’t forget to specify the correct location of mapping_info.txt in the fastq_info entry
   * Don’t forget to specify your mapping directory (should be different from demultiplex directory)
   * Go to the align section and make sure proper reference genome paths are specified
   * Go to the contacts section and make sure proper chrom sizes and cut site files are specified.
   * Note that the _contacts -> call -> call_params_ section cut site parameters need to be in a specific order. The “restriction-sites” parameter can be specified once for each enzyme used (with a cut site location file), but the “restriction-enzymes” parameter must have the same order and must be specified once for each enzyme used.
3. Run [`scripts/A02a_prepare_mapping.sh`](scripts/A02a_prepare_mapping.sh)
   * Don’t forget to specify the correct location of [`txt/mapping_config_snm3C.yml`](txt/mapping_config_snm3C.yml)
   * This is fast
4. Run [`scripts/A02b_run_mapping.sh`](scripts/A02b_run_mapping.sh)
   * This should be submitted with qsub (for SGE, not sure command for other systems)
   * Make sure that path to mapping_scripts.txt is correct (depends on what you named your mapping directory)
   * Don’t forget to make sure that you run a job array of the correct length - (384)x(# of plates)

# snm3C-seq QC analysis

To compute QC statistics for snm3C-seq data using the protocol above, look at the QC statistics file for an individual cell. Use the following equations to compute important QC statistics for this cell.

### Contacts
* cis/trans ratio

  $\frac{(pairs\_enzyme\_intra1kb + pairs\_enzymeless\_intra1kb)}{(pairs\_enzyme\_inter + pairs\_enzymeless\_inter)}$

  Typically, we want this to be greater than 1.

* contacts/pairs ratio

  $\frac{(pairs\_enzyme\_intra1kb + pairs\_enzymeless\_intra1kb + pairs\_enzyme\_inter + pairs\_enzymeless\_inter)}{trimmed\_pairs}$

  We want this to be >0.15

* 3C contact duplicate rate

  $(pairs\_dup\_rate)$

  We want this to be as low as possible. Haven’t benchmarked very stringently, but should be <0.4

### Methylation

* Read duplicates

  $\frac{dupsifter\_removed\_duplicate\_mates}{dupsifter\_input\_mapped\_mates}$

  We want this to be as low as possible. Haven’t benchmarked very stringently, but should be <0.4

* CpG methylation fraction of lambda phage DNA

  $(mCG/CG\_chrL)$

  We want this to be as low as possible, <0.01

* CpH methylation fraction of lambda phage DNA

  $(mCH/CH\_chrL)$

  We want this to be as low as possible, <0.01

* CpG methylation fraction of genomic DNA

  $(mCG/CG\_global)$

  We want this to be >0.5 usually

* CH methylation fraction of genomic DNA

  $(mCH/CH\_global)$

  Usually <0.01 except in neurons

* Presence of non-bisulfite-converted R1 reads (i.e. contamination)

  $\frac{R1\_contam\_fail}{(R1\_contam\_fail + R1\_contam\_pass)}$

  We want this to be as low as possible (<0.01)

* Presence of non-bisulfite-converted R2 reads (i.e. contamination)

  $\frac{R2\_contam\_fail}{(R2\_contam\_fail + R2\_contam\_pass)}$

  We want this to be as low as possible (<0.01)


