#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=25G
#SBATCH --job-name=04_deconvolution_DWLS_MeanRatio_top25
#SBATCH -c 1
#SBATCH -o logs/04_deconvolution_DWLS_MeanRatio_top25.txt
#SBATCH -e logs/04_deconvolution_DWLS_MeanRatio_top25.txt
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
module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 04_deconvolution_DWLS.R MeanRatio_top25 ../../processed-data/08_bulk_deconvolution/markers_MeanRatio_top25.txt

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.1.0
## available from http://research.libd.org/slurmjobs/
