#!/bin/bash

cd /scratch/jmroger/dna-processing-v1/data
gzip ERR1722799_1.fastq
gzip ERR1722799_2.fastq
gzip ERR1722800_1.fastq
gzip ERR1722800_2.fastq
gzip ERR1722803_1.fastq
gzip ERR1722803_2.fastq
gzip ERR1722804_1.fastq
gzip ERR1722804_2.fastq
mv ERR1722799_1.fastq.gz ERR1722799.1.fastq.gz
mv ERR1722799_2.fastq.gz ERR1722799.2.fastq.gz
mv ERR1722800_1.fastq.gz ERR1722800.1.fastq.gz
mv ERR1722800_2.fastq.gz ERR1722800.2.fastq.gz
mv ERR1722803_1.fastq.gz ERR1722803.1.fastq.gz
mv ERR1722803_2.fastq.gz ERR1722803.2.fastq.gz
mv ERR1722804_1.fastq.gz ERR1722804.1.fastq.gz
mv ERR1722804_2.fastq.gz ERR1722804.2.fastq.gz
