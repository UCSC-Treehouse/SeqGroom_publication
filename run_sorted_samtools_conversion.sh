#!/bin/bash

# Jackie Roger
# run_samtools_sort.sh
# Sorts bam file by read name.

# sample name is passed in as an argument
sample=$1
log=/scratch/jmroger/bam-processing/${sample}/sam/${sample}_name-sorted_samtools_output.log
og_bam=/private/groups/treehouse/staging-unrestricted/SeqGroom/primary/original/${sample}/${sample}.bam

# run samtools sort
/private/home/jmroger/SeqGroom/bin/samtools sort \
-n \
-o /scratch/jmroger/bam-processing/${sample}/sam/${sample}_name-sorted.bam \
${og_bam} \
2>&1 | tee ${log}

# get time at start
time_start=`date +"%T"`

# run samtools fastq conversion
/private/home/jmroger/SeqGroom/bin/samtools fastq \
-1 /scratch/jmroger/bam-processing/${sample}/sam/${sample}_name-sorted_1.fastq \
-2 /scratch/jmroger/bam-processing/${sample}/sam/${sample}_name-sorted_2.fastq \
-s /dev/null \
/scratch/jmroger/bam-processing/${sample}/sam/${sample}_name-sorted.bam \
2>&1 | tee ${log}

# get time at finish
time_finish=`date +"%T"`

# write times to output log
echo "Start time: " ${time_start} >> ${log}
echo "Finish time: " ${time_finish} >> ${log}