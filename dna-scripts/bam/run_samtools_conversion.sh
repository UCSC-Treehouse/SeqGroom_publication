#!/bin/bash

# Set up
sample=$1
bam=${sample}.bam
R1=${sample}.ungroomed.R1.fastq
R2=${sample}.ungroomed.R2.fastq
wd=/scratch/jmroger/dna-processing-v1/data/bam/${sample}
samtools=/scratch/jmroger/dna-processing-v1/bin/samtools-1.9/samtools
log=${wd}/logs/${sample}-ungroomed-conversion.log

# Run samtools conversion
cd ${wd}/ungroomed
$samtools fastq \
-1 $R1 \
-2 $R2 \
$bam \
2>&1 | tee $log

gzip $R1
gzip $R2
