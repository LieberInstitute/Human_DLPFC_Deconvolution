#!/bin/bash
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -N GTEx_deconvolution
#$ -o logs/01_GTEx_deconvolution.txt
#$ -e logs/01_GTEx_deconvolution.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R/4.3

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_GTEx_deconvolution.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/