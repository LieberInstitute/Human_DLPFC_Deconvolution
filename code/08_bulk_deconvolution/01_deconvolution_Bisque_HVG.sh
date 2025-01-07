#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=25G
#SBATCH --job-name=01_deconvolution_Bisque_HVG
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --mail-type=ALL
#SBATCH --array=1-10%20

## Define loops and appropriately subset each variable for the array task ID
all_HVG=(10 20 30 40 50 60 70 80 90 100)
HVG=${all_HVG[$(( $SLURM_ARRAY_TASK_ID / 1 % 10 ))]}

## Explicitly pipe script output to a log
log_path=logs/01_deconvolution_Bisque_HVG${HVG}_${SLURM_ARRAY_TASK_ID}.txt

{
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
module load conda_R/4.4

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_deconvolution_Bisque_HVG.R HVG${HVG} ../../processed-data/06_marker_genes/09_HVGs/HVG${HVG}.txt

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/

