#!/bin/bash

# Set up
filetype=$1
sample=$2
groomstatus=$3
R1=$4
R2=$5
wd=/scratch/jmroger/processing-v2/data/${filetype}/${sample}
STAR=/scratch/jmroger/processing-v2/bin/STAR-2.7.1a/bin/Linux_x86_64/STAR
log=${wd}/logs/${sample}-${groomstatus}-STAR.log

# Run STAR
cd ${wd}/$groomstatus
$STAR --runMode alignReads \
--outSAMtype BAM Unsorted \
--readFilesCommand zcat \
--genomeDir /scratch/jmroger/processing-v2/references/star-genome-index \
--outFileNamePrefix ${sample}- \
--readFilesIn ${wd}/${groomstatus}/$R1 ${wd}/${groomstatus}/$R2 \
--quantMode TranscriptomeSAM \
--outSAMattributes NH HI AS NM MD \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--sjdbScore 1 \
2>&1 | tee $log
