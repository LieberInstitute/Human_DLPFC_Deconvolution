#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=2G,h_vmem=2G,h_fsize=400G
#$ -N update_permissions_human_dlpfc_deconvolution
#$ -o logs/update_permissions.txt
#$ -e logs/update_permissions.txt
#$ -m e
#$ -hold_jid move_bulk_round1

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


MAINDIR="/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution"

## Set new permissions
find ${MAINDIR} -type d -exec nfs4_setfacl -a "A:g:hickslab@cm.cluster:RWX" {} \;
find ${MAINDIR} -type d -exec nfs4_setfacl -a "A:gfdi:hickslab@cm.cluster:RWX" {} \;
find ${MAINDIR} -type f -exec nfs4_setfacl -a "A:g:hickslab@cm.cluster:RW" {} \;

find ${MAINDIR} -type d -exec nfs4_setfacl -a "A:g:lieber_lcolladotor@cm.cluster:RWX" {} \;
find ${MAINDIR} -type d -exec nfs4_setfacl -a "A:gfdi:lieber_lcolladotor@cm.cluster:RWX" {} \;
find ${MAINDIR} -type f -exec nfs4_setfacl -a "A:g:lieber_lcolladotor@cm.cluster:RW" {} \;

find ${MAINDIR} -type d -exec nfs4_setfacl -a "A:g:lieber_marmaypag@cm.cluster:RWX" {} \;
find ${MAINDIR} -type d -exec nfs4_setfacl -a "A:gfdi:lieber_marmaypag@cm.cluster:RWX" {} \;
find ${MAINDIR} -type f -exec nfs4_setfacl -a "A:g:lieber_marmaypag@cm.cluster:RW" {} \;

# find ${MAINDIR} -type d -exec nfs4_setfacl -a "A:g:lieber@cm.cluster:RX" {} \;
# find ${MAINDIR} -type d -exec nfs4_setfacl -a "A:gfdi:lieber@cm.cluster:RX" {} \;
# find ${MAINDIR} -type f -exec nfs4_setfacl -a "A:g:lieber@cm.cluster:R" {} \;

## To move away from lieber_jaffe
chgrp lieber_lcolladotor -R ${MAINDIR}

## For setting the group sticky bit
find ${MAINDIR} -type d | xargs chmod g+s

## Check settings
nfs4_getfacl ${MAINDIR}

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/





