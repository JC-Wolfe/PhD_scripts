#!/bin/bash
# File: BG3_Proximal.sh
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -N BG3_Proximal
#$ -l mem_free=10G
#$ -m be
#$ -M jw18713@essex.ac.uk
#$ -o output_BG3_Proximal.txt

./BG3_Expression_Histograms_Proximal.R
