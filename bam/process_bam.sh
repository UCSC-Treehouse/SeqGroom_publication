#!/bin/bash

# Set up
sample=$1
wd=/scratch/jmroger/processing-v2/data/bam/$sample
scripts=/scratch/jmroger/processing-v2/scripts
log=${wd}/logs/${sample}-SUMMARY.log
cd $wd
mkdir logs ungroomed sorted groomed
echo `date` "Bam processing summary for $sample" >> $log
echo >> $log
u=${sample}.bam
cp $u ungroomed
cp $u sorted
mv $u groomed
cd ungroomed
bash ${scripts}/get_bam_info.sh bam $sample ungroomed original $u || true
echo >> $log

# Samtools conversion
cd ${wd}/ungroomed
echo `date` "Samtools conversion" >> $log
u1=${sample}.ungroomed.R1.fastq.gz
u2=${sample}.ungroomed.R2.fastq.gz
bash ${scripts}/bam/run_samtools_conversion.sh $sample || true
bash ${scripts}/get_fastq_info.sh bam $sample ungroomed fastq $u1 $u2 || true
echo >> $log

# Name-sorted samtools conversion
cd ${wd}/sorted
echo `date` "Sorted samtools conversion" >> $log
s=${sample}.sorted.bam
s1=${sample}.sorted.R1.fastq.gz
s2=${sample}.sorted.R2.fastq.gz
bash ${scripts}/bam/run_sorted_samtools_conversion.sh $sample || true
bash ${scripts}/get_fastq_info.sh bam $sample sorted fastq $s1 $s2 || true
echo >> $log

# Seqgroom conversion
cd ${wd}/groomed
echo `date` "SeqGroom conversion" >> $log
docker run --rm -v ${wd}/groomed/:/data \
-e file1=$u \
jackieroger/seqgroom \
2>&1 | tee ${wd}/logs/${sample}-SeqGroom-conversion.log
g1=${sample}.SeqGroomed.R1.fastq.gz
g2=${sample}.SeqGroomed.R2.fastq.gz
bash ${scripts}/get_fastq_info.sh bam $sample groomed fastq $g1 $g2 || true
echo >> $log

# Run rna-seq analysis software on ungroomed fastq files
bash ${scripts}/run_tools.sh bam $sample ungroomed $u1 $u2 || true

# Run rna-seq analysis software on sorted fastq files
bash ${scripts}/run_tools.sh bam $sample sorted $s1 $s2 || true

# Run rna-seq analysis software on groomed fastq files
bash ${scripts}/run_tools.sh bam $sample groomed $g1 $g2 || true

# Finishing up
echo `date` "Deleting extra file copies and compressing subdirectories" >> $log
cd ${wd}/ungroomed
mv $u $wd
rm ${wd}/sorted/${u} ${wd}/groomed/${u}
cd $wd
tar -zcvf ungroomed.tar.gz ungroomed
tar -zcvf sorted.tar.gz sorted
tar -zcvf groomed.tar.gz groomed
rm -rf ungroomed sorted groomed
echo `date` "Processing done!" >> $log
