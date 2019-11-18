#!/bin/bash

# Set up
filetype=$1
sample=$2
groomstatus=$3
info=$4
vcf=$5
wd=/scratch/jmroger/processing-v2/data/${filetype}/$sample
log=${wd}/logs/${sample}-SUMMARY.log

# Get info
cd ${wd}/$groomstatus
vcf_size=`du -sh $vcf`
echo `date` "$groomstatus $info vcf size: $vcf_size" >> $log
vcf_md=`md5sum $bam`
echo `date` "$groomstatus $info vcf md5sum: $vcf_md" >> $log
