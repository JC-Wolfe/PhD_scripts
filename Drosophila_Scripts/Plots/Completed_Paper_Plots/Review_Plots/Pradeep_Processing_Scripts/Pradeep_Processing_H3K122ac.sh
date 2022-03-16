#!/bin/bash
# File: merlinwrapper.sh
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -N k122bowtie
#$ -M jw18713@essex.ac.uk
#$ -m be

# H3K122

# Input
bowtie2 -p 12 -x hg38 -1 Input_5_EKDL210001409-1a-34_HWYY5DSXY_L3_1.fq \
-2 Input_5_EKDL210001409-1a-34_HWYY5DSXY_L3_2.fq \
-S H3K122_input.sam

samtools view -S -b H3K122_input.sam > H3K122_input.bam

# abc_new
bowtie2 -p 12 -x hg38 -1 K122Ac_Abc_New_EKDL210001409-1a-31_HWYY5DSXY_L3_1.fq \
-2 K122Ac_Abc_New_EKDL210001409-1a-31_HWYY5DSXY_L3_2.fq \
-S H3K122_abc_new.sam

samtools view -S -b H3K122_abc_new.sam > H3K122_abc_new.bam

mkdir macs2

macs2 callpeak -t H3K122_abc_new.bam \
-c ../input/H3K122_input.bam \
-f BAM -g 3049315783 \
-n H3K122_abc_new -B \
--outdir macs2 2> macs2/H3K122_abc_new.log

# ac_rs
bowtie2 -p 12 -x hg38 -1 K122Ac_RS_EKDL210001409-1a-29_HWYY5DSXY_L3_1.fq \
-2 K122Ac_RS_EKDL210001409-1a-29_HWYY5DSXY_L3_2.fq \
-S H3K122_rs.sam

samtools view -S -b H3K122_rs.sam > H3K122_rs.bam

mkdir macs2

macs2 callpeak -t H3K122_rs.bam \
-c ../input/H3K122_input.bam \
-f BAM -g 3049315783 \
-n H3K122_rs -B \
--outdir macs2 2> macs2/H3K122_rs.log

#ac_ab_old
bowtie2 -p 12 -x hg38 -1 K122Ac_ab_old_EKDL210001409-1a-33_HWYY5DSXY_L3_1.fq \
-2 K122Ac_ab_old_EKDL210001409-1a-33_HWYY5DSXY_L3_2.fq \
-S H3K122_ab_old.sam

samtools view -S -b H3K122_ab_old.sam > H3K122_ab_old.bam

mkdir macs2

macs2 callpeak -t H3K122_ab_old.bam \
-c ../input/H3K122_input.bam \
-f BAM -g 3049315783 \
-n H3K122_ab_old -B \
--outdir macs2 2> macs2/H3K122_ab_old.log

# Combined
mkdir combined

bowtie2 -p 12 -x hg38 -1 k122ac_abc_new/K122Ac_Abc_New_EKDL210001409-1a-31_HWYY5DSXY_L3_1.fq,\
k122ac_ab_old/K122Ac_ab_old_EKDL210001409-1a-33_HWYY5DSXY_L3_1.fq \
-2 k122ac_abc_new/K122Ac_Abc_New_EKDL210001409-1a-31_HWYY5DSXY_L3_2.fq,\
k122ac_ab_old/K122Ac_ab_old_EKDL210001409-1a-33_HWYY5DSXY_L3_2.fq \
-S combined/H3K122ac_combined.sam

cd combined

samtools view -S -b H3K122ac_combined.sam > H3K122ac_combined.bam

mkdir macs2

macs2 callpeak -t H3K122ac_combined.bam \
-c ../input/H3K122_input.bam \
-f BAM -g 3049315783 \
-n H3K122ac_combined -B \
--outdir macs2 2> macs2/H3K122ac_combined.log
