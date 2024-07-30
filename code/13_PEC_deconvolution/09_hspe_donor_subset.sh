#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=8G
#SBATCH --job-name=09_hspe_donor_subset
#SBATCH -c 1
#SBATCH -t 4:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-998%20

## Define loops and appropriately subset each variable for the array task ID
all_n_donors=(3 4 5 7 10 14 19 27 37 52)
n_donors=${all_n_donors[$(( $SLURM_ARRAY_TASK_ID / 100 % 10 ))]}

all_run_num=($(seq 1 100))
run_num=${all_run_num[$(( $SLURM_ARRAY_TASK_ID / 1 % 100 ))]}

## Explicitly pipe script output to a log
log_path=logs/09_hspe_donor_subset_${n_donors}donors_run${run_num}_${SLURM_ARRAY_TASK_ID}.txt

{
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

## Edit with your job command
Rscript 09_hspe_donor_subset.R --n_donors ${n_donors} --run_num ${run_num}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/

