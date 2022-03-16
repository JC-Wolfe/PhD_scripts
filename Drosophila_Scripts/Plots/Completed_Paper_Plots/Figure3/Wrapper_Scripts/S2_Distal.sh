#!/bin/bash
# File: S2_Distal.sh
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -N S2_Distal
#$ -l mem_free=10G
#$ -m be
#$ -M jw18713@essex.ac.uk
#$ -o output_S2_Distal.txt

./S2_Expression_Histograms_Distal.R
