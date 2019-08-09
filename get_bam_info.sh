#!/bin/bash

# Set up
filetype=$1
sample=$2
groomstatus=$3
info=$4
bam=$5
samtools=/scratch/jmroger/processing-v2/bin/samtools-1.9/samtools
wd=/scratch/jmroger/processing-v2/data/${filetype}/$sample
log=${wd}/logs/${sample}-SUMMARY.log

# Get info
cd ${wd}/$groomstatus
bam_size=`du -sh $bam`
echo `date` "$groomstatus $info bam size: $bam_size" >> $log
bam_md=`md5sum $bam`
echo `date` "$groomstatus $info bam md5sum: $bam_md" >> $log
header=`$samtools view -H $bam | grep -v '@SQ'`
echo `date` "$groomstatus $info bam header info:" >> $log
echo `date` "$header" >> $log
