#!/bin/bash -l
#SBATCH --output=logs/06_DE_library-combo.txt
#SBATCH --error=logs/06_DE_library-combo.txt
#SBATCH --partition=shared
#SBATCH --job-name=DE_library-combo
# 
#SBATCH --mem=75GB

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.3

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 06_DE_library-combo.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
