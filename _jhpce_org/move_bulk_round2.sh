#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=2G,h_vmem=2G,h_fsize=400G
#$ -N move_bulk_round1
#$ -o logs/move_bulk_round1.$TASK_ID.txt
#$ -e logs/move_bulk_round1.$TASK_ID.txt
#$ -m e
#$ -t 1-2
#$ -tc 2

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## List current modules for reproducibility
module list

## Locate directory
ORIGINALDIR=$(awk "NR==${SGE_TASK_ID}" bulk_round1_dirs.txt)
echo "Processing directory ${ORIGINALDIR}"
date

BASEDIR=$(basename ${ORIGINALDIR})
ORIGINALHOME=$(dirname ${ORIGINALDIR})
NEWDIR="/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution/raw-data/bulkRNA/${BASEDIR}"

## Determine amount of data to transfer
du -sk --apparent-size ${ORIGINALDIR}/ | awk '{$1=$1/(1024^3); print $1, "TB";}'

## List owners of all the files in the original path
find ${ORIGINALDIR}/ -exec ls -l {} \; | grep -v total | tr -s ' ' | cut -d ' ' -f3,9-

## Copy from dcl01 to dcs04
rsync -rltgvh --chown=:lieber_lcolladotor ${ORIGINALDIR}/ ${NEWDIR}/

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/