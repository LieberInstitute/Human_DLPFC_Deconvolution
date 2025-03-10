#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=40G
#SBATCH --job-name=08_bisque_donor_subset
#SBATCH -c 4
#SBATCH -t 3:00:00
#SBATCH -o logs/08_bisque_donor_subset_%a.txt
#SBATCH -e /dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution/code/13_PEC_deconvolution/logs/08_bisque_donor_subset_%a.txt
#SBATCH --array=2-10%10

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

## Edit with your job command
Rscript 08_bisque_donor_subset.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.2
## available from http://research.libd.org/slurmjobs/
