#!/bin/bash

# Set up
filetype=$1
sample=$2
groomstatus=$3
R1=$4
R2=$5
wd=/scratch/jmroger/dna-processing-v1/data/${filetype}/${sample}
scripts=/scratch/jmroger/dna-processing-v1/scripts
log=${wd}/logs/${sample}-SUMMARY.log
cd ${wd}/$groomstatus

# Run cutadapt on fastqs
echo `date` "Running cutadapt on $groomstatus fastqs" >> $log
R1_c=${sample}-cutadapt.R1.fastq.gz
R2_c=${sample}-cutadapt.R2.fastq.gz
bash ${scripts}/run_cutadapt.sh $filetype $sample $groomstatus $R1 $R2 $R1_c $R2_c || true
bash ${scripts}/get_fastq_info.sh $filetype $sample $groomstatus cutadapt $R1_c $R2_c || true
echo >> $log

# Run bwa
echo `date` "Running bwa on $groomstatus fastqs" >> $log
bash ${scripts}/run_bwa.sh $filetype $sample $groomstatus $R1 $R2 || true
R_b=${sample}.bwa.bam
bash ${scripts}/get_bam_info.sh $filetype $sample $groomstatus bwa $R_b || true
echo >> $log

# Run deepvariant
echo `date` "Running deepvariant on $groomstatus fastqs" >> $log
bash ${scripts}/run_deepvariant.sh $filetype $sample $groomstatus $R1 $R2 || true
R_d=${sample}.deepvariant.vcf
bash ${scripts}/get_vcf_info.sh $filetype $sample $groomstatus deepvariant $R_d || true
echo >> $log
