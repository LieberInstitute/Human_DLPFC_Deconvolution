#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --job-name=07_deconvolution_CIBERSORTx_HVG
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --mail-type=ALL
#SBATCH --array=1,4-9%10

## Define loops and appropriately subset each variable for the array task ID
all_HVG=(10 20 30 40 50 60 70 80 90 100)
HVG=${all_HVG[$(( $SLURM_ARRAY_TASK_ID / 1 % 10 ))]}

## Explicitly pipe script output to a log
log_path=logs/07_deconvolution_CIBERSORTx_HVG${HVG}_${SLURM_ARRAY_TASK_ID}.txt

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

## Load the module
module load cibersortx

## List current modules for reproducibility
module list

# Define input and output directories
in_dir=/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution/processed-data/08_bulk_deconvolution/07_deconvolution_CIBERSORTx_prep
out_dir=/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution/processed-data/08_bulk_deconvolution/07_deconvolution_CIBERSORTx_prep/output_HVG${HVG}

## make directory
# mkdir ${out_dir}

ref=${in_dir}/DLPFC_sc_counts-HVG${HVG}.txt
mix=${in_dir}/DLPFC_bulk_counts.txt

# Run CIBERSORTxFractions
singularity exec \
    -B ${in_dir}:/src/data \
    -B ${out_dir}:/src/outdir \
    /jhpce/shared/libd/core/cibersortx/04_04_2020/fractions_latest.sif /src/CIBERSORTxFractions \
        --username louise.huuki@libd.org \
        --token $CIBERSORT_TOKEN \
        --single_cell TRUE \
        --refsample $ref \
        --mixture $mix \
        --fraction 0 \
        --rmbatchSmode TRUE \
        --outdir ${out_dir}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/

