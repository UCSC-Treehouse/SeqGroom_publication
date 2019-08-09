#!/bin/bash

STAR=/scratch/jmroger/processing-v2/bin/STAR-2.7.1a/bin/Linux_x86_64/STAR
wd=/scratch/jmroger/processing-v2/references

cd $wd

$STAR --runMode genomeGenerate \
--genomeDir ${wd}/star-genome-index \
--genomeFastaFiles ${wd}/GRCh38.fa \
--sjdbGTFfile ${wd}/gencode.v31.annotation.gtf
