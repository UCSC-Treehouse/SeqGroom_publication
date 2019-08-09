#!/bin/bash

RSEM=/scratch/jmroger/processing-v2/bin/RSEM-1.3.1
wd=/scratch/jmroger/processing-v2/references/rsem-genome-index
refs_wd=/scratch/jmroger/processing-v2/references

cd $wd

${RSEM}/rsem-prepare-reference \
--gtf ${refs_wd}/gencode.v31.annotation.gtf \
${refs_wd}/GRCh38.fa \
hg38
