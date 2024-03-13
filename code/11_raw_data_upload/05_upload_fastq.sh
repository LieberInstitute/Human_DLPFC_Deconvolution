#!/bin/bash

#SBATCH -p transfer
#SBATCH -c 1
#SBATCH --mem=3G
#SBATCH -t 3-00:00:00
#SBATCH --job-name=05_upload_fastq
#SBATCH -o ../../processed-data/11_raw_data_upload/05_upload_fastq.log
#SBATCH -e ../../processed-data/11_raw_data_upload/05_upload_fastq.log
#SBATCH --open-mode=append

set -e

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"

fastq_dir=$(git rev-parse --show-toplevel)/raw-data/bulkRNA/flat_dir
ascp \
    -i ~/aspera.openssh \
    -QT \
    -l100m \
    -k1 \
    -d \
    $fastq_dir \
    $ASP_DEST_DIR # privately sourced environment variable for destination dir

echo "**** Job ends ****"
date
