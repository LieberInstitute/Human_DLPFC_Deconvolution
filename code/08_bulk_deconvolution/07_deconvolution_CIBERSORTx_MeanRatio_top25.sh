#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=25G
#SBATCH --job-name=07_deconvolution_CIBERSORTx_MeanRatio_top25
#SBATCH -c 1
#SBATCH -o logs/07_deconvolution_CIBERSORTx_MeanRatio_top25.txt
#SBATCH -e logs/07_deconvolution_CIBERSORTx_MeanRatio_top25.txt
#SBATCH --mail-type=ALL

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the module
module load cibersortx

## List current modules for reproducibility
module list

# Define input and output directories
in_dir=/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution/processed-data/08_bulk_deconvolution/07_deconvolution_CIBERSORTx_prep
out_dir=/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution/processed-data/08_bulk_deconvolution/07_deconvolution_CIBERSORTx_prep/output_MeanRatio_top25

# Run CIBERSORTxFractions
singularity exec \
-B ${in_dir}:/src/data \
-B ${out_dir}/:/src/outdir \
/jhpce/shared/libd/core/cibersortx/04_04_2020/fractions_latest.sif /src/CIBERSORTxFractions \
--username louise.huuki@libd.org \
--token $CIBERSORT_TOKEN \
--single_cell TRUE \
--refsample ${in_dir}/DLPFC_sc_counts-MeanRatio_top25.txt \
--mixture ${in_dir}/DLPFC_bulk_counts.txt \
--fraction 0 \
--rmbatchSmode TRUE \
--outdir ${out_dir}


echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.1.0
## available from http://research.libd.org/slurmjobs/
