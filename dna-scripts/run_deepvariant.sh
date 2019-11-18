#!/bin/bash

# Set up
filetype=$1
sample=$2
groomstatus=$3
R1=$4
R2=$5
wd=/scratch/jmroger/dna-processing-v1/data/${filetype}/${sample}
log=${wd}/logs/${sample}-${groomstatus}-deepvariant.log

# Run deepvariant
cd ${wd}/$groomstatus
INPUT_DIR=/scratch/jmroger/dna-processing-v1
OUTPUT_DIR=/scratch/jmroger/dna-processing-v1
docker run \
-v "${INPUT_DIR}":"/input" \
-v "${OUTPUT_DIR}:/output" \
gcr.io/deepvariant-docker/deepvariant:0.8.0 \
/opt/deepvariant/bin/run_deepvariant \
--model_type=WGS \
--ref=/input/references/hg19.fa \
--reads=/input/data/${filetype}/${sample}/${groomstatus}/${sample}.bwa.bam \
--sample_name=${sample} \
--output_vcf=/output/data/${filetype}/${sample}/${groomstatus}/${sample}.deepvariant.vcf.gz \
2>&1 | tee $log
