#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH --job-name=02_sra_biosample
#SBATCH -o ../../processed-data/11_raw_data_upload/02_sra_biosample.log
#SBATCH -e ../../processed-data/11_raw_data_upload/02_sra_biosample.log

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
Rscript 02_sra_biosample.R

echo "**** Job ends ****"
date
