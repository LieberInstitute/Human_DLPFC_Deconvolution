#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 10:00
#SBATCH --job-name=04_prepare_fastqs
#SBATCH -o ../../processed-data/11_raw_data_upload/04_prepare_fastqs.log
#SBATCH -e ../../processed-data/11_raw_data_upload/04_prepare_fastqs.log

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.3
Rscript 04_prepare_fastqs.R

echo "**** Job ends ****"
date
