#!/bin/bash -l
#SBATCH --output=logs/05_DE_library-prep.txt
#SBATCH --error=logs/05_DE_library-prep.txt
#SBATCH --partition=shared
#SBATCH --job-name=05_DE_library-prep
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
Rscript 05_DE_library-prep.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/