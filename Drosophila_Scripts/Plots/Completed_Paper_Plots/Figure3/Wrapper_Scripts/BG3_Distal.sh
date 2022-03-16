#!/bin/bash
# File: BG3_Distal.sh
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -N BG3_Distal
#$ -l mem_free=10G
#$ -m be
#$ -M jw18713@essex.ac.uk
#$ -o output_BG3_Distal.txt

./BG3_Expression_Histograms_Distal.R
