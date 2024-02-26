#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH --job-name=01_deconvolution_Bisque_random_subset
#SBATCH -c 4
#SBATCH -o logs/01_deconvolution_Bisque_random_subset.txt
#SBATCH -e logs/01_deconvolution_Bisque_random_subset.txt

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
module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_deconvolution_Bisque_random_subset.R

echo "**** Job ends ****"
date
