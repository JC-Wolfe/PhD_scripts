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

# H4K16ac September

# Experimental Files
bowtie2 -p 12 -x hg38 -1 T1_H4K16ac_03_EKDL200003000-1a-N703_HGNGGDSXY_L3_1.fq \
-2 T1_H4K16ac_03_EKDL200003000-1a-N703_HGNGGDSXY_L3_2.fq \
-S H4K16ac_sep.sam

# Input
bowtie2 -p 12 -x hg38 -1 T1_Rb_IgG_10_EKDL200003000-1a-N710_HGNGGDSXY_L3_1.fq \
-2 T1_Rb_IgG_10_EKDL200003000-1a-N710_HGNGGDSXY_L3_2.fq \
-S Rb_IgG_sep.sam

samtools view -S -b Rb_IgG_sep.sam > Rb_IgG_sep.bam

samtools view -S -b H4K16ac_sep.sam > H4K16ac_sep.bam

mkdir macs2

macs2 callpeak -t H4K16ac_sep.bam \
-c Rb_IgG_sep.bam \
-f BAM -g 3049315783 \
-n H4K16ac_sep -B \
--outdir macs2 2> macs2/H4K16ac_sep.log

# H4K16ac August

# Experimental Files
bowtie2 -p 12 -x hg38 -1 H4K16ac_EKDL200002317-1a-N704_HCYLJDSXY_L3_1.fq \
-2 H4K16ac_EKDL200002317-1a-N704_HCYLJDSXY_L3_2.fq \
-S H4K16ac_aug.sam

# Input
bowtie2 -p 12 -x hg38 -1 IgG_EKDL200002317-1a-N708_HCYLJDSXY_L3_1.fq \
-2 IgG_EKDL200002317-1a-N708_HCYLJDSXY_L3_2.fq \
-S Rb_IgG_aug.sam


samtools view -S -b Rb_IgG_aug.sam > Rb_IgG_aug.bam

samtools view -S -b H4K16ac_aug.sam > H4K16ac_aug.bam

mkdir macs2

macs2 callpeak -t H4K16ac_aug.bam \
-c Rb_IgG_aug.bam \
-f BAM -g 3049315783 \
-n H4K16ac_aug -B \
--outdir macs2 2> macs2/H4K16ac_aug.log

#H4K16 combined
# Experimental Files
bowtie2 -p 12 -x hg38 -1 aug_2020/H4K16ac_EKDL200002317-1a-N704_HCYLJDSXY_L3_1.fq,sep_2020/T1_H4K16ac_03_EKDL200003000-1a-N703_HGNGGDSXY_L3_1.fq \
-2 aug_2020/H4K16ac_EKDL200002317-1a-N704_HCYLJDSXY_L3_2.fq,sep_2020/T1_H4K16ac_03_EKDL200003000-1a-N703_HGNGGDSXY_L3_2.fq \
-S combined/H4K16ac_combined.sam

# Inputs
bowtie2 -p 12 -x hg38 -1 aug_2020/IgG_EKDL200002317-1a-N708_HCYLJDSXY_L3_1.fq,sep_2020/T1_Rb_IgG_10_EKDL200003000-1a-N710_HGNGGDSXY_L3_1.fq \
-2 aug_2020/IgG_EKDL200002317-1a-N708_HCYLJDSXY_L3_2.fq,sep_2020/T1_Rb_IgG_10_EKDL200003000-1a-N710_HGNGGDSXY_L3_2.fq \
-S combined/H4K16ac_input_combined.sam

cd combined

samtools view -S -b H4K16ac_combined.sam > H4K16ac_combined.bam

samtools view -S -b H4K16ac_input_combined.sam > H4K16ac_input_combined.bam

macs2 callpeak -t H4K16ac_combined.bam \
-c H4K16ac_input_combined.bam \
-f BAM -g 3049315783 \
-n H4K16ac_combined -B \
--outdir macs2 2> macs2/H4K16ac_combined.log

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

macs2 callpeak -t 3K122_ab_old.bam \
-c ../input/H3K122_input.bam \
-f BAM -g 3049315783 \
-n 3K122_ab_old -B \
--outdir macs2 2> macs2/3K122_ab_old.log
