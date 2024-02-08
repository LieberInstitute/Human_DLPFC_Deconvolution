#!/bin/bash
<<<<<<< HEAD
#$ -cwd
#$ -l mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -N 01_deconvolution_Bisque_FULL
#$ -o logs/01_deconvolution_Bisque_FULL.txt
#$ -e logs/01_deconvolution_Bisque_FULL.txt
#$ -m e
=======
#SBATCH -p shared
#SBATCH --mem=25G
#SBATCH --job-name=01_deconvolution_Bisque_FULL
#SBATCH -c 1
#SBATCH -o logs/01_deconvolution_Bisque_FULL.txt
#SBATCH -e logs/01_deconvolution_Bisque_FULL.txt
#SBATCH --mail-type=ALL

set -e
>>>>>>> 6a6608cb57c6d7a265a8624f14ea7cb4bcdf3e02

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

<<<<<<< HEAD
## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.3.x
=======
## Load the R module
module load conda_R/4.3
>>>>>>> 6a6608cb57c6d7a265a8624f14ea7cb4bcdf3e02

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_deconvolution_Bisque.R FULL

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.1.0
## available from http://research.libd.org/slurmjobs/
