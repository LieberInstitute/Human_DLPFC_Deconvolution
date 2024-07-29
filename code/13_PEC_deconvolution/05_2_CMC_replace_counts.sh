#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=60G
#SBATCH --job-name=05_2_CMC_replace_counts
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o logs/05_2_CMC_replace_counts.txt
#SBATCH -e logs/05_2_CMC_replace_counts.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R/4.3.x

## List current modules for reproducibility
module list

Rscript 05_2_CMC_replace_counts.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.2
## available from http://research.libd.org/slurmjobs/
