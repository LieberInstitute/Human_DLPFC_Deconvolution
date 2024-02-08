#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=25G
#SBATCH --job-name=05_deconvolution_hspe_1vALL_top25
#SBATCH -c 1
#SBATCH -o logs/05_deconvolution_hspe_1vALL_top25.txt
#SBATCH -e logs/05_deconvolution_hspe_1vALL_top25.txt
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
module load conda_R/4.3

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 05_deconvolution_hspe.R 1vALL_top25

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.1.0
## available from http://research.libd.org/slurmjobs/
