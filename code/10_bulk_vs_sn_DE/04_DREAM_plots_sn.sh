#!/bin/bash
#SBATCH -p shared
#SBATCH --mem-per-cpu=5G
#SBATCH --job-name=04_DREAM_plots_sn
#SBATCH -c 1
#SBATCH -o logs/04_DREAM_plots_sn.txt
#SBATCH -e logs/04_DREAM_plots_sn.txt
#SBATCH --mail-type=ALL

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 04_DREAM_plots_sn.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 0.99.0
## available from http://research.libd.org/slurmjobs/