#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=5G
#SBATCH --job-name=02_deconvolution_MuSiC_cell_size
#SBATCH -c 1
#SBATCH -o logs/02_deconvolution_MuSiC_cell_size_%a.txt
#SBATCH -e logs/02_deconvolution_MuSiC_cell_size_%a.txt
#SBATCH --array=1-3%3

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

Rscript 02_deconvolution_MuSiC_cell_size.R

echo "**** Job ends ****"
date
