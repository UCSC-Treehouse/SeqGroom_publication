#!/bin/bash

# Jackie Roger
# bam_processing_data_mining.sh
# Mines relevant data from bam file conversion processing.


# ----------------------------------------------------------------------
# SAMPLE INFO
# ----------------------------------------------------------------------

# sample name is passed in as an argument
sample=$1
log=${sample}.log
og_bam=/private/groups/treehouse/staging-unrestricted/SeqGroom/primary/original/${sample}/${sample}.bam

echo "Sample name: ${sample}" >> ${log}

echo -e "\n" >> ${log}



########################################################################
# BAM INFO (og, sorted)
########################################################################

# ----------------------------------------------------------------------
# ORIGINAL BAM INFO
# ----------------------------------------------------------------------

echo "ORIGINAL BAM INFO" >> ${log}

echo -e "\n" >> ${log}

# find the size of the original bam file
og_bam_size=`du -h ${og_bam} | cut -f1`
echo "Size of original bam file: ${og_bam_size}" >> ${log}

echo -e "\n" >> ${log}

# find all lines starting with @HD
HD_lines=`/private/home/jmroger/SeqGroom/bin/samtools view -H ${og_bam} | grep ^@HD`
echo "@HD line(s):" >> ${log}
echo ${HD_lines} >> ${log}

echo -e "\n" >> ${log}

# find  all lines starting with @RG
HD_lines=`/private/home/jmroger/SeqGroom/bin/samtools view -H ${og_bam} | grep ^@RG`
echo "@RG line(s):" >> ${log}
echo ${RG_lines} >> ${log}

echo -e "\n" >> ${log}

# find all lines starting with @PG
HD_lines=`/private/home/jmroger/SeqGroom/bin/samtools view -H ${og_bam} | grep ^@PG`
echo "@PG line(s):" >> ${log}
echo ${PG_lines} >> ${log}

echo -e "\n" >> ${log}

# find all lines starting with @CO
HD_lines=`/private/home/jmroger/SeqGroom/bin/samtools view -H ${og_bam} | grep ^@CO`
echo "@CO line(s):" >> ${log}
echo ${CO_lines} >> ${log}

echo -e "\n\n" >> ${log}

# ----------------------------------------------------------------------
# NAME-SORTED BAM INFO
# ----------------------------------------------------------------------

echo "NAME-SORTED BAM INFO" >> ${log}

echo -e "\n" >> ${log}

# name-sorted bam file
sorted_bam="/scratch/jmroger/bam-processing/${sample}/sam/${sample}.bam"

# find the size of the original bam file
sorted_bam_size=`du -h ${sorted_bam} | cut -f1`
echo "Size of name-sorted bam file: ${sorted_bam_size}" >> ${log}

echo -e "\n\n" >> ${log}



########################################################################
# CONVERSION INFO (samtools, sorted samtools, seqgroom)
########################################################################

# ----------------------------------------------------------------------
# SAMTOOLS CONVERSION INFO
# ----------------------------------------------------------------------

echo "SAMTOOLS CONVERSION INFO" >> ${log}

echo -e "\n" >> ${log}

# samtools fastq files
samtools_fastq_1="scratch/jmroger/bam-processing/${sample}/sam/${sample}_1.fastq"
samtools_fastq_2="scratch/jmroger/bam-processing/${sample}/sam/${sample}_2.fastq"

# samtools log
samtools_log="/scratch/jmroger/bam-processing/${sample}/sam/${sample}_samtools_output.log"

# start and finish times
samtools_start_time=`grep "Start time" ${samtools_log}`
samtools_finish_time=`grep "Finish time" ${samtools_log}`
echo ${samtools_start_time} >> ${log}
echo ${samtools_finish_time} >> ${log}

# size and line count of fastqs
samtools_fastq_1_size=`gzip ${samtools_fastq_1} | du -h`
samtools_fastq_2_size=`gzip ${samtools_fastq_2} | du -h`
samtools_fastq_1_linecount=`wc -l ${samtools_fastq_1}`
samtools_fastq_2_linecount=`wc -l ${samtools_fastq_2}`
echo "R1 fastq file size: ${samtools_fastq_1_size}" >> ${log}
echo "R2 fastq file size: ${samtools_fastq_2_size}" >> ${log}
echo "R1 fastq file linecount: ${samtools_fastq_1_linecount}" >> ${log}
echo "R2 fastq file linecount: ${samtools_fastq_2_linecount}" >> ${log}

echo -e "\n\n" >> ${log}

# ----------------------------------------------------------------------
# NAME-SORTED SAMTOOLS CONVERSION INFO
# ----------------------------------------------------------------------

echo "NAME-SORTED SAMTOOLS CONVERSION INFO" >> ${log}

echo -e "\n" >> ${log}

# sorted samtools fastq files
samtools_fastq_1="scratch/jmroger/bam-processing/${sample}/sam/${sample}_name-sorted_1.fastq"
samtools_fastq_2="scratch/jmroger/bam-processing/${sample}/sam/${sample}_name-sorted_2.fastq"

# samtools log
samtools_log="/scratch/jmroger/bam-processing/${sample}/sam/${sample}_name-sorted_samtools_output.log"

# start and finish times
samtools_start_time=`grep "Start time" ${samtools_log}`
samtools_finish_time=`grep "Finish time" ${samtools_log}`
echo ${samtools_start_time} >> ${log}
echo ${samtools_finish_time} >> ${log}

# size and line count of fastqs
samtools_fastq_1_size=`gzip ${samtools_fastq_1} | du -h`
samtools_fastq_2_size=`gzip ${samtools_fastq_2} | du -h`
samtools_fastq_1_linecount=`wc -l ${samtools_fastq_1}`
samtools_fastq_2_linecount=`wc -l ${samtools_fastq_2}`
echo "R1 fastq file size: ${samtools_fastq_1_size}" >> ${log}
echo "R2 fastq file size: ${samtools_fastq_2_size}" >> ${log}
echo "R1 fastq file linecount: ${samtools_fastq_1_linecount}" >> ${log}
echo "R2 fastq file linecount: ${samtools_fastq_2_linecount}" >> ${log}

echo -e "\n\n" >> ${log}

# ----------------------------------------------------------------------
# SEQGROOM CONVERSION INFO
# ----------------------------------------------------------------------

echo "SEQGROOM CONVERSION INFO" >> ${log}

echo -e "\n" >> ${log}

# seqgroom fastq files
seqgroom_fastq_1="scratch/jmroger/bam-processing/${sample}/seq/${sample}.btfv9.R1.fastq.gz"
seqgroom_fastq_2="scratch/jmroger/bam-processing/${sample}/seq/${sample}.btfv9.R2.fastq.gz"

# samtools log
seqgroom_log="/scratch/jmroger/bam-processing/${sample}/seq/${sample}_seqgroom_output.log"

# start and finish times
seqgroom_start_time=`grep "Start time" ${seqgroom_log}`
seqgroom_finish_time=`grep "Finish time" ${seqgroom_log}`
echo ${seqgroom_start_time} >> ${log}
echo ${seqgroom_finish_time} >> ${log}

# size and line count of fastqs
seqgroom_fastq_1_size=`du -h ${seqgroom_fastq_1}`
seqgroom_fastq_2_size=`du -h ${seqgroom_fastq_2}`
seqgroom_fastq_1_linecount=`gunzip ${seqgroom_fastq_1} | wc -l`
seqgroom_fastq_2_linecount=`gunzip ${seqgroom_fastq_2} | wc -l`
echo "R1 fastq file size: ${samtools_fastq_1_size}" >> ${log}
echo "R2 fastq file size: ${samtools_fastq_2_size}" >> ${log}
echo "R1 fastq file linecount: ${samtools_fastq_1_linecount}" >> ${log}
echo "R2 fastq file linecount: ${samtools_fastq_2_linecount}" >> ${log}

echo -e "\n\n" >> ${log}



########################################################################
# CGL FAILURE INFO (samtools, sorted samtools, seqgroom)
########################################################################

# ----------------------------------------------------------------------
# SAMTOOLS CGL FAILURE INFO
# ----------------------------------------------------------------------

echo "SAMTOOLS CGL FAILURE INFO" >> ${log}

echo -e "\n" >> ${log}

sam_cgl_log`=${sample}_samtools_cgl_output.log`

sam_cgl_error=`grep ^ERROR ${sam_cgl_log}`
echo "Error message:" >> ${log}
echo ${sam_cgl_error} >> ${log}
echo -e "\n" >> ${log}

sam_cgl_warning=`grep ^WARNING ${sam_cgl_log}`
echo "Warning message:" >> ${log}
echo ${sam_cgl_warning} >> ${log}

echo -e "\n\n" >> ${log}

# ----------------------------------------------------------------------
# NAME-SORTED SAMTOOLS CGL FAILURE INFO
# ----------------------------------------------------------------------

echo "NAME-SORTED SAMTOOLS CGL FAILURE INFO" >> ${log}

echo -e "\n" >> ${log}

sorted_sam_cgl_log`=${sample}_name-sorted_samtools_cgl_output.log`

sorted_sam_cgl_error=`grep ^ERROR ${sorted_sam_cgl_log}`
echo "Error message:" >> ${log}
echo ${sorted_sam_cgl_error} >> ${log}
echo -e "\n" >> ${log}

sorted_sam_cgl_warning=`grep ^WARNING ${sorted_sam_cgl_log}`
echo "Warning message:" >> ${log}
echo ${sorted_sam_cgl_warning} >> ${log}

echo -e "\n\n" >> ${log}

# ----------------------------------------------------------------------
# SEQGROOM CGL FAILURE INFO
# ----------------------------------------------------------------------

echo "SEQGROOM CGL BAM INFO" >> ${log}

echo -e "\n" >> ${log}

seq_cgl_log`=${sample}_seqgroom_cgl_output.log`

seq_cgl_error=`grep ^ERROR ${seq_cgl_log}`
echo "Error message:" >> ${log}
echo ${seq_cgl_error} >> ${log}
echo -e "\n" >> ${log}

seq_cgl_warning=`grep ^WARNING ${seq_cgl_log}`
echo "Warning message:" >> ${log}
echo ${seq_cgl_warning} >> ${log}

echo -e "\n\n" >> ${log}



########################################################################
# CGL BAM INFO (samtools, sorted samtools, seqgroom)
########################################################################

# ----------------------------------------------------------------------
# SAMTOOLS CGL BAM INFO
# ----------------------------------------------------------------------

echo "SAMTOOLS CGL BAM INFO" >> ${log}

echo -e "\n" >> ${log}

# samtools cgl bam file
sam_cgl_bam="/scratch/jmroger/bam-processing/${sample}/sam/${sample}_samtools_cgl.bam"

# find the size of the samtools cgl file
sam_cgl_bam_size=`du -h ${sam_cgl_bam} | cut -f1`
echo "Size of samtools cgl bam file: ${sam_cgl_bam_size}" >> ${log}

echo -e "\n\n" >> ${log}

# ----------------------------------------------------------------------
# NAME-SORTED SAMTOOLS CGL BAM INFO
# ----------------------------------------------------------------------

echo "NAME-SORTED SAMTOOLS CGL BAM INFO" >> ${log}

echo -e "\n" >> ${log}

# name-sorted samtools cgl bam file
sorted_sam_cgl_bam="/scratch/jmroger/bam-processing/${sample}/sam/${sample}_name-sorted_samtools_cgl.bam"

# find the size of the name-sortedsamtools cgl file
sorted_sam_cgl_bam_size=`du -h ${sorted_sam_cgl_bam} | cut -f1`
echo "Size of name-sorted samtools cgl bam file: ${sorted_sam_cgl_bam_size}" >> ${log}

echo -e "\n\n" >> ${log}

# ----------------------------------------------------------------------
# SEQGROOM CGL BAM INFO
# ----------------------------------------------------------------------

echo "SEQGROOM CGL BAM INFO" >> ${log}

echo -e "\n" >> ${log}

# seqgroom cgl bam file
seq_cgl_bam="/scratch/jmroger/bam-processing/${sample}/seq/${sample}_seqgroom_cgl.bam"

# find the size of the seqgroom cgl file
seq_cgl_bam_size=`du -h ${sam_cgl_bam} | cut -f1`
echo "Size of seqgroom cgl bam file: ${seq_cgl_bam_size}" >> ${log}

echo -e "\n\n" >> ${log}



########################################################################
# CGL BAM & RSEM INFO (samtools, sorted samtools, seqgroom)
########################################################################

# ----------------------------------------------------------------------
# SAMTOOLS CGL BAM & RSEM INFO
# ----------------------------------------------------------------------

echo "SAMTOOLS CGL BAM & RSEM INFO" >> ${log}

echo -e "\n" >> ${log}

# samtools cgl bam file
sam_cgl_bam="/scratch/jmroger/bam-processing/${sample}/sam/${sample}_samtools_cgl.bam"

# find the size of the samtools cgl file
sam_cgl_bam_size=`du -h ${sam_cgl_bam} | cut -f1`
echo "Size of samtools cgl bam file: ${sam_cgl_bam_size}" >> ${log}

echo -e "\n\n" >> ${log}

# ----------------------------------------------------------------------
# NAME-SORTED SAMTOOLS CGL BAM & RSEM INFO
# ----------------------------------------------------------------------

echo "NAME-SORTED SAMTOOLS CGL & RSEM INFO" >> ${log}

echo -e "\n" >> ${log}

# name-sorted samtools cgl bam file
sorted_sam_cgl_bam="/scratch/jmroger/bam-processing/${sample}/sam/${sample}_name-sorted_samtools_cgl.bam"

# find the size of the name-sortedsamtools cgl file
sorted_sam_cgl_bam_size=`du -h ${sorted_sam_cgl_bam} | cut -f1`
echo "Size of name-sorted samtools cgl bam file: ${sorted_sam_cgl_bam_size}" >> ${log}

echo -e "\n\n" >> ${log}

# ----------------------------------------------------------------------
# SEQGROOM CGL BAM & RSEM INFO
# ----------------------------------------------------------------------

echo "SEQGROOM CGL BAM & RSEM INFO" >> ${log}

echo -e "\n" >> ${log}

# seqgroom cgl bam file
seq_cgl_bam="/scratch/jmroger/bam-processing/${sample}/seq/${sample}_seqgroom_cgl.bam"

# find the size of the seqgroom cgl file
seq_cgl_bam_size=`du -h ${sam_cgl_bam} | cut -f1`
echo "Size of seqgroom cgl bam file: ${seq_cgl_bam_size}" >> ${log}

echo -e "\n\n" >> ${log}

