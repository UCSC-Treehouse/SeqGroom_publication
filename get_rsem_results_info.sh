#!/bin/bash

# Set up
filetype=$1
sample=$2
groomstatus=$3
info=$4
rsem_results=$5
wd=/scratch/jmroger/processing-v2/data/${filetype}/$sample
log=${wd}/logs/${sample}-SUMMARY.log

# Get info
cd ${wd}/$groomstatus
rsem_size=`du -sh $rsem_results`
echo `date` "$groomstatus $info results size: $rsem_size" >> $log
rsem_md=`md5sum $rsem_results`
echo `date` "$groomstatus $info results md5sum: $rsem_md" >> $log
