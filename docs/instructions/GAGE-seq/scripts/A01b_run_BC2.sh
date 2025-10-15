#!/bin/bash
#$ -cwd
#$ -o logs/A01b_run_BC2.$JOB_ID.$TASK_ID
#$ -j y
#$ -N A01b_run_BC2
#$ -l h_data=10G,h_rt=20:00:00
#$ -pe shared 3
#$ -t 1-1:1

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

source ~/.bashrc

conda activate map3C_snakemake # <-

SNAKE="demultiplex/results/demultiplex_BC2_scripts.txt"

ID=$SGE_TASK_ID

RUN=`head -${ID} $SNAKE | tail -1`

time bash $RUN

echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `date `
echo " "
