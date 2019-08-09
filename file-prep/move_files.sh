#!/bin/bash

s=$1

og=/scratch/jmroger/freely-available-files/fastq
new=/scratch/jmroger/processing-v2/fastq

mv ${og}/${s}/${s}_1.fastq.gz ${new}/${s}/${s}.R1.fastq.gz
mv ${og}/${s}/${s}_2.fastq.gz ${new}/${s}/${s}.R2.fastq.gz
