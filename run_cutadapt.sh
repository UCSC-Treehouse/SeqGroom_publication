#!/bin/bash

# Set up
filetype=$1
sample=$2
groomstatus=$3
R1=$4
R2=$5
R1_c=$6
R2_c=$7
wd=/scratch/jmroger/processing-v2/data/${filetype}/${sample}
cutadapt=/scratch/jmroger/processing-v2/bin/cutadapt
log=${wd}/logs/${sample}-${groomstatus}-cutadapt.log

# Run cutadapt
cd ${wd}/$groomstatus
$cutadapt \
-a AGATCGGAAGAG -m 35 -A AGATCGGAAGAG \
-o $R1_c -p $R2_c $R1 $R2 \
2>&1 | tee $log

