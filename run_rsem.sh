#!/bin/bash

# Set up
filetype=$1
sample=$2
groomstatus=$3
bam=$4
wd=/scratch/jmroger/processing-v2/data/${filetype}/${sample}
RSEM=/scratch/jmroger/processing-v2/bin/RSEM-1.3.1
log=${wd}/logs/${sample}-${groomstatus}-RSEM.log

# Run RSEM
cd ${wd}/$groomstatus
${RSEM}/rsem-calculate-expression \
--quiet \
--no-qualities \
--forward-prob 0.5 \
--seed-length 25 \
--fragment-length-mean -1.0 \
--bam \
--paired-end \
$bam \
/scratch/jmroger/processing-v2/references/rsem-genome-index/hg38 \
${sample} \
2>&1 | tee $log
