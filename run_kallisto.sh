#!/bin/bash

# Set up
filetype=$1
sample=$2
groomstatus=$3
R1=$4
R2=$5
wd=/scratch/jmroger/processing-v2/data/${filetype}/${sample}
kallisto=/scratch/jmroger/processing-v2/bin/kallisto/build/src/kallisto
log=${wd}/logs/${sample}-${groomstatus}-kallisto.log

# Run kallisto
cd ${wd}/$groomstatus
$kallisto quant \
-i /scratch/jmroger/processing-v2/references/kallisto-genome-index/transcriptome.idx \
-o ${sample}-kallisto-results \
-b 100 \
--fusion \
$R1 $R2 \
2>&1 | tee $log
