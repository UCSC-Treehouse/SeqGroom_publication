#!/bin/bash

# Jackie Roger
# run_samtools_conversion.sh
# Converts bam to fastq using samtools fastq.

# sample name is passed in as an argument
sample=$1
log=/scratch/jmroger/bam-processing/${sample}/sam/${sample}_samtools_output.log
og_bam=/private/groups/treehouse/staging-unrestricted/SeqGroom/primary/original/${sample}/${sample}.bam

# get time at start
time_start=`date +"%T"`

# run samtools fastq conversion
/private/home/jmroger/SeqGroom/bin/samtools fastq \
-1 /scratch/jmroger/bam-processing/${sample}/sam/${sample}_1.fastq \
-2 /scratch/jmroger/bam-processing/${sample}/sam/${sample}_2.fastq \
${og_bam} \
2>&1 | tee ${log}

# get time at finish
time_finish=`date +"%T"`

# write times to output log
echo "Start time: " ${time_start} >> ${log}
echo "Finish time: " ${time_finish} >> ${log}