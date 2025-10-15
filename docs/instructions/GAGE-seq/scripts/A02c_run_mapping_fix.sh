#!/bin/bash
#$ -cwd
#$ -o logs/A02c_run_mapping_fix.$JOB_ID.$TASK_ID
#$ -j y
#$ -N A02c_run_mapping_fix
#$ -l h_data=10G,h_rt=20:00:00
#$ -pe shared 2
#$ -t 1-7:1

ulimit -c 0
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

source ~/.bashrc

conda activate map3C_snakemake # <-

SNAKE="txt/mapping_cmd_fix.txt"

ID=$SGE_TASK_ID

RUN=`head -${ID} $SNAKE | tail -1`

time bash $RUN

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "
