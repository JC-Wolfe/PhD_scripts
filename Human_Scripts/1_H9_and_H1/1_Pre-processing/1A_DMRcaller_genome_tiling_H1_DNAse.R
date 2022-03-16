# Importing libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("genomation")
#BiocManager::install("DMRcaller")

library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(rtracklayer)
library(genomation)
library(DMRcaller)

# Chain for lift overs
chain <- import.chain("/home/jw18713/Human_Datasets/hg19ToHg38.over.chain")

# Getting the hg38 genome
genome <- getSeq(BSgenome.Hsapiens.UCSC.hg38)

# Getting my own functions
source("~/Genome_construction_functions.R")

# First I need to build a 100bp object to store my results in
# This will be my master object to add scores to going forward
H1_tiled <- empty_tiled_genome(Hsapiens, binsize = 100, chr_range = seq(1,23))

# Next I need an empty tiled genome set at every 10 base pairs to average the
# reads over
H1_10bin <- empty_tiled_genome(Hsapiens, binsize = 10, chr_range = seq(1,23))

# Now I need to populate the mcols with track information
setwd("~/Human_Datasets/H1/ChIP/")
folders <- dir()

for (i in seq_along(folders)){
  setwd(paste0("~/Human_Datasets/H1/ChIP/", folders[i]))
  for (j in seq(1,length(dir()))){
    mark <- unlist(strsplit(getwd(), '/'))[length(unlist(strsplit(getwd(), '/')))]
    mark <- paste0(mark, "_", j)
    print(paste0(mark, ": ", i, " of ", length(folders)))
    track <- import(dir()[j])
    mcols(H1_tiled)[length(mcols(H1_tiled)) + 1] <- track_calculator(H1_10bin,
      track)
    names(mcols(H1_tiled))[length(mcols(H1_tiled))] <- mark
  }
}

# Saving the histone modifications and clearing memory
save(H1_tiled, file = "~/Human_Datasets/H1/DMR_caller_method_pt1.Rda")

load("~/Human_Datasets/H1/DMR_caller_method_pt1.Rda")


# Adding accessibility data with ATAC-seq
setwd("~/Human_Datasets/H1/ATAC-seq")

# Loading in ATAC-seq tracks, extraCols due to MACS2 peak calling
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
ATAC1 <- import("GSM2264818_H1_0_1.filterBL.bed", format = "BED", extraCols = extraCols_narrowPeak)
ATAC2 <- import("GSM2264819_H1_0_2.filterBL.bed", format = "BED", extraCols = extraCols_narrowPeak)

# Lifting over from hg19 to hg39
ATAC1 <- unlist(liftOver(ATAC1, chain))
ATAC2 <- unlist(liftOver(ATAC2, chain))

# Running DMRcaller based track_calculator to average reads in 10bp bins over
# 100bp windows
ATAC1_score <- track_calculator(H1_10bin, ATAC1)
ATAC2_score <- track_calculator(H1_10bin, ATAC2)
ATAC1_sv <- track_calculator(H1_10bin, ATAC1, target = "signalValue")
ATAC2_sv <- track_calculator(H1_10bin, ATAC2, target = "signalValue")

H1_tiled$ATAC_scores <- apply(cbind(ATAC1_score, ATAC2_score), 1, FUN = mean)
H1_tiled$ATAC_signalValues <- apply(cbind(ATAC1_sv, ATAC2_sv), 1, FUN = mean)

# Saving file with ATAC-seq completed
save(H1_tiled, file = "DMR_caller_method_after_ATAC_seq.Rda")

# Bisulfite sequencing
setwd("~/Human_Datasets/H1/Bisulfite")

folders <- dir()
for (i in seq_along(folders)) {
  setwd(paste0("~/Human_Datasets/H1/Bisulfite/",folders[i]))
  name <- folders[i]
  t1 <- import(dir()[1])
  t2 <- import(dir()[2])
  t1_score <- Bisulfite_calculator(H1_tiled, t1)
  t2_score <- Bisulfite_calculator(H1_tiled, t2)
  mcols(H1_tiled)[length(mcols(H1_tiled)) + 1] <- apply(cbind(t1_score, t2_score), 1, FUN = mean)
  names(mcols(H1_tiled))[length(mcols(H1_tiled))] <- name
}

# Saving Bisulfite data
save(H1_tiled, file = "~/Human_Datasets/H1/DMRcaller_tiled_noSTARR.Rda")

# Adding DNAse data
load("~/Human_Datasets/H1/DMRcaller_tiled_noSTARR.Rda")
DNAse <- import("~/Human_Datasets/H1/DNAse/ENCFF405MMH.bigWig")
H1_tiled$DNAse <- track_calculator(H1_10bin, DNAse)
save(H1_tiled, file = "~/Human_Datasets/H1/H1_after_DNAse.Rda")

# Normalising the columns by min of max normalisation
for (c in seq_along(mcols(H1_tiled))){
  max_value = max(mcols(H1_tiled)[c][,1])
  min_value = min(mcols(H1_tiled)[c][,1])
  mcols(H1_tiled)[c][,1] <- (mcols(H1_tiled)[c][,1] - min_value) / (max_value - min_value)
}

save(H1_tiled, file = "~/Human_Datasets/H1/DMRcaller_norm_tiled.Rda")
