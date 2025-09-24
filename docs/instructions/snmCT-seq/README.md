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
STAR  --runMode genomeGenerate \
--runThreadN 20 \
--genomeDir /path/to/output/index/ \
--genomeFastaFiles /path/to/ref.fa \
--sjdbGTFfile /path/to/gtf \
 --sjdbOverhang 149 \
--limitGenomeGenerateRAM 50000000000
```

Also, you will need to download the chromosome size files for your reference genome, as well as a GTF for gene locations


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
2. Update [`txt/demultiplex_config_snmCT.yml`](txt/demultiplex_config_snmCT.yml) with:
   * The path to [`txt/plate_info.txt`](txt/plate_info.txt) goes in the fastq_info entry
   * The path to demultiplex directory - this will include your demultiplexed FASTQ files and no mapping results - goes in the output_directory entry
   * snmCTseq barcodes FASTA (the proper file is provided at [`txt/random_index_v2.multiplex.fa`](txt/random_index_v2.multiplex.fa))
3. Run [`scripts/A01a_prepare_demultiplex.sh`](scripts/A01a_prepare_demultiplex.sh)
   * Make sure to update the path to [`txt/demultiplex_config.yml`](txt/demultiplex_config_snmCT.yml)
   * This is fast
4. Run [`scripts/A01b_run_demultiplex.sh`](scripts/A01b_run_demultiplex.sh)
   * This should be submitted with qsub (for SGE, not sure command for other systems)
   * Don’t forget to make sure that you run a job array of the correct length (number of plates) - i.e. if you have 24 plates, the qsub parameter should be -t 1-24:1
   * A01b should finish in a couple hours for each plate. Check to make sure that there is a file for each plate that has this format: `demultiplex/results/{plate}/{plate}_demultiplex_stats.txt`

# Running map3C (mapping)

1. Update [`txt/mapping_info.txt`](txt/mapping_info.txt)
   * TSV where first column is well name (underscores are allowed), second column is the whole path to R1 FASTQ, and third column is the whole path to R2 FASTQ
   * Should have (384)x(# of plates) lines
2. Update [`txt/mapping_config_snmCT.yml`](txt/mapping_config_snmCT.yml)
   * Don’t forget to specify the correct location of mapping_info.txt in the fastq_info entry
   * Don’t forget to specify your mapping directory (should be different from demultiplex directory)
   * Go to the _dna -> align -> biscuit -> reference_path_ section and make sure proper reference genome path is specified
   * Go to the _dna -> postprocess -> call_contacts_params_ section and make sure proper chrom sizes path is specified
   * Go to the _rna -> align -> star -> reference_path section_ and make sure proper reference genome path is specified
   * Go to the _rna -> postprocess -> featurecounts_pe_params_ and featurecounts_pe_params sections and make sure proper GTF path is specified
3. Run [`scripts/A02a_prepare_mapping.sh`](scripts/A02a_prepare_mapping.sh)
   * Don’t forget to specify the correct location of [`txt/mapping_config_snmCT.yml`](txt/mapping_config_snmCT.yml)
   * This is fast
4. Run [`scripts/A02b_run_mapping.sh`](scripts/A02b_run_mapping.sh)
   * This should be submitted with qsub (for SGE, not sure command for other systems)
   * Make sure that path to mapping_scripts.txt is correct (depends on what you named your mapping directory)
   * Don’t forget to make sure that you run a job array of the correct length - (384)x(# of plates)
