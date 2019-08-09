#!/bin/bash

# Set up
filetype=$1
sample=$2
groomstatus=$3
R1=$4
R2=$5
wd=/scratch/jmroger/processing-v2/data/${filetype}/${sample}
scripts=/scratch/jmroger/processing-v2/scripts
log=${wd}/logs/${sample}-SUMMARY.log
cd ${wd}/$groomstatus

# Run cutadapt on fastqs
echo `date` "Running cutadapt on $groomstatus fastqs" >> $log
R1_c=${sample}-cutadapt.R1.fastq.gz
R2_c=${sample}-cutadapt.R2.fastq.gz
bash ${scripts}/run_cutadapt.sh $filetype $sample $groomstatus $R1 $R2 $R1_c $R2_c || true
bash ${scripts}/get_fastq_info.sh $filetype $sample $groomstatus cutadapt $R1_c $R2_c || true
echo >> $log

# Run STAR & RSEM on fastqs
echo `date` "Running STAR on $groomstatus fastqs" >> $log
bash ${scripts}/run_star.sh $filetype $sample $groomstatus $R1 $R2 || true
R_s=${sample}-Aligned.toTranscriptome.out.bam
bash ${scripts}/get_bam_info.sh $filetype $sample $groomstatus star $R_s || true
echo `date` "Running RSEM on STAR-aligned $groomstatus reads" >> $log
bash ${scripts}/run_rsem.sh $filetype $sample $groomstatus $R_s || true
R_r=${sample}.genes.results
bash ${scripts}/get_rsem_results_info.sh $filetype $sample $groomstatus rsem $R_r || true
echo >> $log

# Run kallisto on fastqs
echo `date` "Running kallisto on $groomstatus fastqs" >> $log
bash ${scripts}/run_kallisto.sh $filetype $sample $groomstatus $R1 $R2 || true
R_k=${sample}-kallisto-results/abundance.tsv
bash ${scripts}/get_kallisto_results_info.sh $filetype $sample $groomstatus kallisto $R_k || true
echo >> $log
