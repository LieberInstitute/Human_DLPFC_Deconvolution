#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -N spatial_size_QC
#$ -o logs/02_spatial_size_QC.txt
#$ -e logs/02_spatial_size_QC.txt
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
Rscript 02_spatial_size_QC.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
