#!/bin/bash

# Set up
sample=$1
cd /scratch/jmroger/dna-processing-v1/data/fastq/
mkdir $sample
mv ${sample}.R1.fastq.gz $sample
mv ${sample}.R2.fastq.gz $sample
wd=/scratch/jmroger/dna-processing-v1/data/fastq/$sample
scripts=/scratch/jmroger/dna-processing-v1/scripts
log=${wd}/logs/${sample}-SUMMARY.log
cd $wd
mkdir logs ungroomed groomed
echo `date` "Fastq processing summary for $sample" >> $log
echo >> $log
u1=${sample}.R1.fastq.gz
u2=${sample}.R2.fastq.gz
cp $u1 $u2 ungroomed
mv $u1 $u2 groomed
cd ungroomed
bash ${scripts}/get_fastq_info.sh fastq $sample ungroomed fastq $u1 $u2 || true
echo >> $log

# Groom fastq files using seqgroom
cd ${wd}/groomed
echo `date` "Running SeqGroom on ungroomed fastqs" >> $log
docker run --rm -v ${wd}/groomed/:/data \
-e file1=$u1 \
-e file2=$u2 \
jackieroger/seqgroom \
2>&1 | tee ${wd}/logs/${sample}-SeqGroom.log
g1=${sample}.SeqGroomed.R1.fastq.gz
g2=${sample}.SeqGroomed.R2.fastq.gz
bash ${scripts}/get_fastq_info.sh fastq $sample groomed fastq $g1 $g2 || true
echo >> $log

# Run tools on ungroomed fastq files
bash ${scripts}/run_tools.sh fastq $sample ungroomed $u1 $u2 || true

# Run tools on groomed fastq files
bash ${scripts}/run_tools.sh fastq $sample groomed $g1 $g2 || true

# Finishing up
echo `date` "Deleting extra file copies and compressing subdirectories" >> $log
rm $u1 $u2
cd ${wd}/ungroomed
mv $u1 $u2 $wd
cd $wd
tar -zcvf ungroomed.tar.gz ungroomed
tar -zcvf groomed.tar.gz groomed
rm -rf ungroomed groomed
echo `date` "Processing done!" >> $log
