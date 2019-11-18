#!/bin/bash

# Set up
filetype=$1
sample=$2
groomstatus=$3
R1=$4
R2=$5
wd=/scratch/jmroger/dna-processing-v1/data/${filetype}/${sample}
bwa=/scratch/jmroger/dna-processing-v1/bin/bwa/bwa
samtools=/scratch/jmroger/dna-processing-v1/bin/samtools-1.9/samtools
log=${wd}/logs/${sample}-${groomstatus}-bwa.log

# Run bwa
cd ${wd}/$groomstatus
$bwa mem /scratch/jmroger/dna-processing-v1/references/hg19.fa \
$R1 $R2 > ${sample}.bwa.sam

# Convert sam to bam and index bam
$samtools view -S -b -h ${sample}.bwa.sam > ${sample}.bwa.bam
$samtools index ${sample}.bwa.bam
