
IgG_Paul
Paul_1_1000
Paul_1_250
Paul_1_500

# H3K4me2
# IgG_Paul
# Input
bowtie2 -p 12 -x hg38 -1 IgG_P_hurd_EKDL200012967-1a-6_HNMHWDSXY_L3_1.fq \
-2 IgG_P_hurd_EKDL200012967-1a-6_HNMHWDSXY_L3_2.fq \
-S H3K4me2_input.sam

samtools view -S -b H3K4me2_input.sam > H3K4me2_input.bam

# Paul_1_1000
bowtie2 -p 12 -x hg38 -1 Paul_1_1000_EKDL200012968-1a-N710_HNMHWDSXY_L2_1.fq \
-2 Paul_1_1000_EKDL200012968-1a-N710_HNMHWDSXY_L2_2.fq \
-S H3K4me2_1000.sam

samtools view -S -b H3K4me2_1000.sam > H3K4me2_1000.bam

mkdir macs2

macs2 callpeak -t H3K4me2_1000.bam \
-c ../IgG_Paul/H3K4me2_input.bam \
-f BAM -g 3049315783 \
-n H3K4me2_1000 -B \
--outdir macs2 2> macs2/H3K4me2_1000.log

# Paul_1_250
bowtie2 -p 12 -x hg38 -1 Paul_1_250_EKDL200012968-1a-N705_HNMHWDSXY_L2_1.fq \
-2 Paul_1_250_EKDL200012968-1a-N705_HNMHWDSXY_L2_2.fq \
-S H3K4me2_250.sam

samtools view -S -b H3K4me2_250.sam > H3K4me2_250.bam

mkdir macs2

macs2 callpeak -t H3K4me2_250.bam \
-c ../IgG_Paul/H3K4me2_input.bam \
-f BAM -g 3049315783 \
-n H3K4me2_250 -B \
--outdir macs2 2> macs2/H3K4me2_250.log

# Paul_1_500
bowtie2 -p 12 -x hg38 -1 Paul_1_500_EKDL200012968-1a-N709_HNMHWDSXY_L2_1.fq \
-2 Paul_1_500_EKDL200012968-1a-N709_HNMHWDSXY_L2_2.fq \
-S H3K4me2_500.sam

samtools view -S -b H3K4me2_500.sam > H3K4me2_500.bam

mkdir macs2

macs2 callpeak -t H3K4me2_500.bam \
-c ../IgG_Paul/H3K4me2_input.bam \
-f BAM -g 3049315783 \
-n H3K4me2_500 -B \
--outdir macs2 2> macs2/H3K4me2_500.log

# H3K4me2 combined
mkdir combined

bowtie2 -p 12 -x hg38 -1 Paul_1_1000/Paul_1_1000_EKDL200012968-1a-N710_HNMHWDSXY_L2_1.fq,\
Paul_1_250/Paul_1_250_EKDL200012968-1a-N705_HNMHWDSXY_L2_1.fq,\
Paul_1_500/Paul_1_500_EKDL200012968-1a-N709_HNMHWDSXY_L2_1.fq \
-2 Paul_1_1000/Paul_1_1000_EKDL200012968-1a-N710_HNMHWDSXY_L2_2.fq,\
Paul_1_250/Paul_1_250_EKDL200012968-1a-N705_HNMHWDSXY_L2_2.fq,\
Paul_1_500/Paul_1_500_EKDL200012968-1a-N709_HNMHWDSXY_L2_2.fq \
-S combined/H3K4me2_combined.sam

cd combined

samtools view -S -b H3K4me2_combined.sam > H3K4me2_combined.bam

mkdir macs2

macs2 callpeak -t H3K4me2_combined.bam \
-c ../IgG_Paul/H3K4me2_input.bam \
-f BAM -g 3049315783 \
-n H3K4me2_combined -B \
--outdir macs2 2> macs2/H3K4me2_combined.log
