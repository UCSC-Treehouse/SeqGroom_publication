#!/bin/bash

# Set up
filetype=$1
sample=$2
groomstatus=$3
info=$4
R1=$5
R2=$6
wd=/scratch/jmroger/processing-v2/data/${filetype}/$sample
log=${wd}/logs/${sample}-SUMMARY.log

# Get info
cd ${wd}/$groomstatus
R1_size=`du -sh ${R1}`
R2_size=`du -sh ${R2}`
echo `date` "$groomstatus $info R1 size: $R1_size" >> $log
echo `date` "$groomstatus $info R2 size: $R2_size" >> $log
R1_lc=`gzip -d -c $R1 | wc -l`
R2_lc=`gzip -d -c $R2 | wc -l`
echo `date` "$groomstatus $info R1 linecount: $R1_lc" >> $log
echo `date` "$groomstatus $info R2 linecount: $R2_lc" >> $log
R1_md=`md5sum $R1`
R2_md=`md5sum $R2`
echo `date` "$groomstatus $info R1 md5sum: $R1_md" >> $log
echo `date` "$groomstatus $info R2 md5sum: $R2_md" >> $log
