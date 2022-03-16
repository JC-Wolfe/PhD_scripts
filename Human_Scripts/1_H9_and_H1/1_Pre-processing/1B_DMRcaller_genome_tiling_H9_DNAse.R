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
H9_tiled <- empty_tiled_genome(Hsapiens, binsize = 100, chr_range = seq(1,23))

# Next I need an empty tiled genome set at every 10 base pairs to average the
# reads over
H9_10bin <- empty_tiled_genome(Hsapiens, binsize = 10, chr_range = seq(1,23))

# Now I need to populate the mcols with track information
setwd("~/Human_Datasets/H9/ChIP/")
folders <- dir()

for (i in seq_along(folders)){
  setwd(paste0("~/Human_Datasets/H9/ChIP/", folders[i]))
  for (j in seq(1,length(dir()))){
    mark <- unlist(strsplit(getwd(), '/'))[length(unlist(strsplit(getwd(), '/')))]
    mark <- paste0(mark, "_", j)
    print(paste0(mark, ": ", i, " of ", length(folders)))
    track <- import(dir()[j])
    mcols(H9_tiled)[length(mcols(H9_tiled)) + 1] <- track_calculator(H9_10bin,
      track)
    names(mcols(H9_tiled))[length(mcols(H9_tiled))] <- mark
  }
}

# Saving the histone modifications and clearing memory
save(H9_tiled, file = "~/Human_Datasets/H9/DMR_caller_method_pt1.Rda")
load("~/Human_Datasets/H9/DMR_caller_method_pt1.Rda")


# Adding accessibility data with ATAC-seq
setwd("~/Human_Datasets/H9/ATAC-seq")

# Loading in ATAC-seq tracks, extraCols due to MACS2 peak calling
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
ATAC1 <- import("GSM2264826_H9_0_1.filterBL.bed", format = "BED", extraCols = extraCols_narrowPeak)
ATAC2 <- import("GSM2264827_H9_0_2.filterBL.bed", format = "BED", extraCols = extraCols_narrowPeak)

# Lifting over from hg19 to hg39
ATAC1 <- unlist(liftOver(ATAC1, chain))
ATAC2 <- unlist(liftOver(ATAC2, chain))

# Running DMRcaller based track_calculator to average reads in 10bp bins over
# 100bp windows
ATAC1_score <- track_calculator(H9_10bin, ATAC1)
ATAC2_score <- track_calculator(H9_10bin, ATAC2)
ATAC1_sv <- track_calculator(H9_10bin, ATAC1, target = "signalValue")
ATAC2_sv <- track_calculator(H9_10bin, ATAC2, target = "signalValue")

H9_tiled$ATAC_scores <- apply(cbind(ATAC1_score, ATAC2_score), 1, FUN = mean)
H9_tiled$ATAC_signalValues <- apply(cbind(ATAC1_sv, ATAC2_sv), 1, FUN = mean)

# Saving file with ATAC-seq completed
save(H9_tiled, file = "DMR_caller_method_after_ATAC_seq.Rda")

# Bisulfite sequencing
setwd("~/Human_Datasets/H9/Bisulfite")

folders <- dir()
for (i in seq_along(folders)) {
  setwd(paste0("~/Human_Datasets/H9/Bisulfite/",folders[i]))
  name <- folders[i]
  t1 <- import(dir()[1])
  t2 <- import(dir()[2])
  t1_score <- Bisulfite_calculator(H9_tiled, t1)
  t2_score <- Bisulfite_calculator(H9_tiled, t2)
  mcols(H9_tiled)[length(mcols(H9_tiled)) + 1] <- apply(cbind(t1_score, t2_score), 1, FUN = mean)
  names(mcols(H9_tiled))[length(mcols(H9_tiled))] <- name
}

# Saving Bisulfite data
save(H9_tiled, file = "~/Human_Datasets/H9/DMRcaller_tiled_noSTARR.Rda")

# Adding DNAse data
load("~/Human_Datasets/H9/DMRcaller_tiled_noSTARR.Rda")
DNAse <- import("~/Human_Datasets/H9/DNAse/ENCFF405MMH.bigWig")
H9_tiled$DNAse <- track_calculator(H9_10bin, DNAse)
save(H9_tiled, file = "~/Human_Datasets/H9/H9_after_DNAse.Rda")

# Adding STARR-seq annotations
STARR <- import("/home/jw18713/Human_Datasets/H9/STARR-seq/ChIP-STARR-seq_enhancers_primed.bed")
STARR <- STARR[STARR$score >= 138]

# STARR-seq liftOver
STARR <- unlist(liftOver(STARR, chain))

save(STARR, file = "~/Human_Datasets/H9/STARR-seq/Processed_STARR.Rda")

# Assigning a 1 to STARR-seq overlaps
overlaps <- findOverlaps(H9_tiled,STARR)
scores <- rep(0, length(H9_tiled))
scores[as.vector(queryHits(overlaps))] <- 1
H9_tiled$STARR_seq_binary <- scores

# Saving the lifted over STARR-seq dataset
save(H9_tiled, file = "~/Human_Datasets/H9/DMRcaller_with_STARR_DNAse.Rda")

# Normalising the columns by min of max normalisation
for (c in seq_along(mcols(H9_tiled))){
  max_value = max(mcols(H9_tiled)[c][,1])
  min_value = min(mcols(H9_tiled)[c][,1])
  mcols(H9_tiled)[c][,1] <- (mcols(H9_tiled)[c][,1] - min_value) / (max_value - min_value)
}

save(H9_tiled, file = "~/Human_Datasets/H9/DMRcaller_norm_tiled_DNAse.Rda")
