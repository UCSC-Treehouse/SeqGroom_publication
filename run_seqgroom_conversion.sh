#!/bin/bash

# Jackie Roger
# run_seqgroom_conversion.sh
# Converts bam to fastq using SeqGroom.

# sample name is passed in as an argument
sample=$1
log=/scratch/jmroger/bam-processing/${sample}/seq/${sample}_seqgroom_output.log
og_bam=/private/groups/treehouse/staging-unrestricted/SeqGroom/primary/original/${sample}/${sample}.bam

# put a copy of the bam into the seq subdirectory
# (btfv9 requires that outputs be in same directory as input)
cp ${og_bam} \
/scratch/jmroger/bam-processing/${sample}/seq

# get time at start
time_start=`date +"%T"`

# run SeqGroom conversion
docker run --rm -v /scratch/jmroger/bam-processing/${sample}/seq/:/data \
-e input=${sample}.bam linhvoyo/btfv9 \
2>&1 | tee ${log}

# get time at finish
time_finish=`date +"%T"`

# write times to output log
echo "Start time: " ${time_start} >> ${log}
echo "Finish time: " ${time_finish} >> ${log}