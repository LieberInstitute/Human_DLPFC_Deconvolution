#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=25G
#SBATCH --job-name=07_deconvolution_CIBERSORTx_test2
#SBATCH -c 1
#SBATCH -o logs/07_deconvolution_CIBERSORTx_test2.txt
#SBATCH -e logs/07_deconvolution_CIBERSORTx_test2.txt
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

## Load the R module
module load cibersortx

## List current modules for reproducibility
module list

## Edit with your job command

singularity exec \
-B /dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution/processed-data/08_bulk_deconvolution/07_deconvolution_CIBERSORTx_prep/tutorial_data/Fig2ab-NSCLC_PBMCs/:/src/data \
-B /dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution/processed-data/08_bulk_deconvolution/07_deconvolution_CIBERSORTx_prep/tutorial_data/Fig2ab-NSCLC_PBMCs/output:/src/outdir \
/jhpce/shared/libd/core/cibersortx/04_04_2020/fractions_latest.sif /src/CIBERSORTxFractions \
--username louise.huuki@libd.org \
--token $CIBERSORT_TOKEN \
--single_cell TRUE \
--refsample /dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution/processed-data/08_bulk_deconvolution/07_deconvolution_CIBERSORTx_prep/tutorial_data/Fig2ab-NSCLC_PBMCs/sce_counts_test.txt \
--mixture /dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution/processed-data/08_bulk_deconvolution/07_deconvolution_CIBERSORTx_prep/tutorial_data/Fig2ab-NSCLC_PBMCs/Fig2b-WholeBlood_RNAseq.txt \
--fraction 0 \
--rmbatchSmode TRUE \
--outdir /dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution/processed-data/08_bulk_deconvolution/07_deconvolution_CIBERSORTx_prep/tutorial_data/Fig2ab-NSCLC_PBMCs/output


echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.1.0
## available from http://research.libd.org/slurmjobs/