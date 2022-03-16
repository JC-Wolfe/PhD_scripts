
bowtie2 -p 12 -x dm6 -1 starr_rep1_1.fastq -2 starr_rep1_2.fastq -S starr_rep1.sam

bowtie2 -p 12 -x dm6 -1 starr_rep2_1.fastq -2 starr_rep2_2.fastq -S starr_rep2.sam

bowtie2 -p 12 -x dm6 -1 starr_input_1.fastq -2 starr_input_2.fastq -S starr_input.sam

samtools view -S -b starr_rep1.sam > starr_rep1.bam

samtools view -S -b starr_rep2.sam > starr_rep2.bam

samtools view -S -b starr_input.sam > starr_input.bam

macs2 callpeak -t starr_rep1.bam \
-c starr_input.bam \
-f BAM -g 142573017 \
-n starr_rep1 -B \
--outdir macs2 2> macs2/starr_rep1_macs2.log

macs2 callpeak -t starr_rep2.bam \
-c starr_input.bam \
-f BAM -g 142573017 \
-n starr_rep2 -B \
--outdir macs2 2> macs2/starr_rep2_macs2.log
