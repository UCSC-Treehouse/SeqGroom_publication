#!/bin/bash

# Set up
sample=$1
og_bam=${sample}.bam
sorted_bam=${sample}.sorted.bam
R1=${sample}.sorted.R1.fastq
R2=${sample}.sorted.R2.fastq
wd=/scratch/jmroger/dna-processing-v1/data/bam/${sample}
samtools=/scratch/jmroger/dna-processing-v1/bin/samtools-1.9/samtools
log=${wd}/logs/${sample}-sorted-conversion.log

# Sort bam
cd ${wd}/sorted
$samtools sort -n \
-o $sorted_bam \
$og_bam \
2>&1 | tee $log

# Run samtools conversion
$samtools fastq \
-1 $R1 \
-2 $R2 \
$sorted_bam \
2>&1 | tee $log

gzip $R1
gzip $R2
