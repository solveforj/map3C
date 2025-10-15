#!/bin/bash
#$ -cwd
#$ -o logs/A02b_run_mapping.$JOB_ID.$TASK_ID
#$ -j y
#$ -N A02b_run_mapping
#$ -l h_data=5G,h_rt=2:00:00
#$ -pe shared 2
#$ -t 1-1:1
ulimit -c 0

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

source ~/.bashrc

conda activate map3C_snakemake # <-

SNAKE="mapping/results/mapping_scripts.txt"

ID=$SGE_TASK_ID

RUN=`head -${ID} $SNAKE | tail -1`

time bash $RUN

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "
