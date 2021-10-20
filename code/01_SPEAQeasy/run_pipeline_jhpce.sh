#!/bin/bash
#$ -l mem_free=40G,h_vmem=40G,h_fsize=800G
#$ -o ./SPEAQeasy_output.log
#$ -e ./SPEAQeasy_output.log
#$ -cwd

#  After running 'install_software.sh', this should point to the directory
#  where SPEAQeasy was installed, and not say "$PWD"
ORIG_DIR=/users/lhuuki/SPEAQeasy

module load nextflow
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

nextflow $ORIG_DIR/main.nf \
    --sample "paired" \
    --reference "hg38" \
    --strand "reverse" \
    --input "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution/raw-data/bulkRNA" \
    --output "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/Human_DLPFC_Deconvolution/processed-data/01_SPEAQeasy/round1_2021-10-19" \
    --annotation "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/SPEAQeasy/Annotation" \
    --experiment "Human_DLPFC_Deconvolution" \
    --strand_mode 'accept' \
    -with-report execution_reports/JHPCE_run.html \
    -profile jhpce \
	-resume

#  Produces a report for each sample tracing the pipeline steps
#  performed (can be helpful for debugging).
#
#  Note that the reports are generated from the output log produced in the above
#  section, and so if you rename the log, you must also pass replace the filename
#  in the bash call below.
echo "Generating per-sample logs for debugging..."
bash $ORIG_DIR/scripts/generate_logs.sh $PWD/SPEAQeasy_output.log
