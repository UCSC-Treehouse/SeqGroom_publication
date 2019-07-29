#!/bin/bash

sample=$1
log=/scratch/jmroger/fastq-processing/${sample}/${sample}_fastq_linecounts.log

ungroomed_fastq_R1=/scratch/jmroger/fastq-processing/${sample}/${sample}_R1.fastq.gz
ungroomed_fastq_R2=/scratch/jmroger/fastq-processing/${sample}/${sample}_R2.fastq.gz
groomed_fastq_R1=/scratch/jmroger/fastq-processing/${sample}/${sample}_btfv9.R1.fastq.gz
groomed_fastq_R2=/scratch/jmroger/fastq-processing/${sample}/${sample}_btfv9.R2.fastq.gz

ungroomed_R1_linecount=`gzip -d -c ${ungroomed_fastq_R1} | wc -l`
ungroomed_R2_linecount=`gzip -d -c ${ungroomed_fastq_R2} | wc -l`
groomed_R1_linecount=`gzip -d -c ${groomed_fastq_R1} | wc -l`
groomed_R2_linecount=`gzip -d -c ${groomed_fastq_R2} | wc -l`

echo "Ungroomed R1 linecount: ${ungroomed_R1_linecount}" >> $log
echo "Ungroomed R2 linecount: ${ungroomed_R2_linecount}" >> $log
echo "Groomed R1 linecount: ${groomed_R1_linecount}" >> $log
echo "Groomed R2 linecount: ${groomed_R2_linecount}" >> $log