#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=3G
#SBATCH --job-name=05_deconvolution_hspe_random_subset
#SBATCH -c 1
#SBATCH -o logs/05_deconvolution_hspe_random_subset_%a.txt
#SBATCH -e logs/05_deconvolution_hspe_random_subset_%a.txt
#SBATCH --array=418,429,853,859,863,882,918,924,935,953,958,964,969,983,988,992,996,999%15

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

#   Circumvent temporary issues with /tmp being full on some compute nodes
export TMPDIR=$MYSCRATCH

## List current modules for reproducibility
module list

Rscript 05_deconvolution_hspe_random_subset.R

echo "**** Job ends ****"
date
