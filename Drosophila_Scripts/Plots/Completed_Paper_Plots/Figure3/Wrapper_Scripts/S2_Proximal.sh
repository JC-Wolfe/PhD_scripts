#!/bin/bash
# File: S2_Proximal.sh
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -N S2_Proximal
#$ -l mem_free=10G
#$ -m be
#$ -M jw18713@essex.ac.uk
#$ -o output_S2_Proximal.txt

./S2_Expression_Histograms_Proximal.R
