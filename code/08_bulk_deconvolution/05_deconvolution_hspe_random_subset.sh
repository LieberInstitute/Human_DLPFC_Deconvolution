#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=3G
#SBATCH --job-name=05_deconvolution_hspe_random_subset
#SBATCH -c 1
#SBATCH -o logs/05_deconvolution_hspe_random_subset_%a.txt
#SBATCH -e logs/05_deconvolution_hspe_random_subset_%a.txt
#SBATCH --array=36,37,38,39,40,41,42,43,44,45,46,48,49,50,55,252,254,407,410,411,414,415,416,417,418,420,421,428,429,430,838,841,844,849,851,853,855,857,859,861,863,865,867,870,872,878,882,886,889,909,911,915,918,921,924,931,933,935,937,939,942,944,950,953,956,958,960,961,964,966,968,969,973,975,977,978,981,983,987,988,992,996,999%15

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

Rscript 05_deconvolution_hspe_random_subset.R

echo "**** Job ends ****"
date
