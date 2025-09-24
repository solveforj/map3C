#!/bin/bash
#$ -cwd
#$ -o logs/A02b_run_mapping.$JOB_ID.$TASK_ID
#$ -j y
#$ -N A02b_run_mapping
#$ -l h_data=1G,h_rt=1:00:00
#$ -pe shared 3
#$ -t 1-768:1

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

source ~/.bashrc

conda activate map3C_snakemake # <-

SNAKE="mapping/results/mapping_scripts.txt"

ID=$SGE_TASK_ID

RUN=`head -${ID} $SNAKE | tail -1`

time bash $RUN

sleep 10

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "
