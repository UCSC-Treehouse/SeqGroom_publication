#!/bin/bash

# Set up
filetype=$1
sample=$2
groomstatus=$3
info=$4
kallisto_results=$5
wd=/scratch/jmroger/processing-v2/data/${filetype}/$sample
log=${wd}/logs/${sample}-SUMMARY.log

# Get info
cd ${wd}/$groomstatus
kallisto_size=`du -sh $kallisto_results`
echo `date` "$groomstatus $info results size: $kallisto_size" >> $log
kallisto_md=`md5sum $kallisto_results`
echo `date` "$groomstatus $info results md5sum: $kallisto_md" >> $log
