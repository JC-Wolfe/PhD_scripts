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
HeLa_tiled <- empty_tiled_genome(Hsapiens, binsize = 100, chr_range = seq(1,23))

# Next I need an empty tiled genome set at every 10 base pairs to average the
# reads over
HeLa_10bin <- empty_tiled_genome(Hsapiens, binsize = 10, chr_range = seq(1,23))

# Now I need to populate the mcols with track information
setwd("~/Human_Datasets/HeLa/ChIP/")
folders <- dir()

for (i in seq_along(folders)){
  setwd(paste0("~/Human_Datasets/HeLa/ChIP/", folders[i]))
  for (j in seq(1,length(dir()))){
    mark <- unlist(strsplit(getwd(), '/'))[length(unlist(strsplit(getwd(), '/')))]
    mark <- paste0(mark, "_", j)
    print(paste0(mark, ": ", i, " of ", length(folders)))
    track <- import(dir()[j])
    mcols(HeLa_tiled)[length(mcols(HeLa_tiled)) + 1] <- track_calculator(HeLa_10bin,
      track)
    names(mcols(HeLa_tiled))[length(mcols(HeLa_tiled))] <- mark
  }
}

# Saving the histone modifications and clearing memory
save(HeLa_tiled, file = "~/Human_Datasets/HeLa/DMR_caller_method_pt1.Rda")
load("~/Human_Datasets/HeLa/DMR_caller_method_pt1.Rda")

# Bisulfite sequencing
setwd("~/Human_Datasets/HeLa/Bisulfite")

bigbedread <- function(x){
  t <- read.csv(x, sep = "\t", header = F)
  tgr <- GRanges(seqnames = t[,1], ranges = IRanges(
    start = t[,2], end = t[,3]))
  tgr$score <- t[,5]
  return(tgr)
}

folders <- dir()
for (i in seq_along(folders)) {
  setwd(paste0("~/Human_Datasets/HeLa/Bisulfite/",folders[i]))
  name <- folders[i]
  t1 <- bigbedread(dir()[1])
  t2 <- bigbedread(dir()[2])
  t1_score <- Bisulfite_calculator(HeLa_tiled, t1)
  t2_score <- Bisulfite_calculator(HeLa_tiled, t2)
  mcols(HeLa_tiled)[length(mcols(HeLa_tiled)) + 1] <- apply(cbind(t1_score, t2_score), 1, FUN = mean)
  names(mcols(HeLa_tiled))[length(mcols(HeLa_tiled))] <- name
}

# Saving Bisulfite data
save(HeLa_tiled, file = "~/Human_Datasets/HeLa/DMRcaller_tiled_with_bisulfite.Rda")

# Adding DNAse data
load("~/Human_Datasets/HeLa/DMRcaller_tiled_noSTARR.Rda")
DNAse <- import("~/Human_Datasets/HeLa/DNAse/ENCFF172PKC.bigWig")
HeLa_tiled$DNAse <- track_calculator(HeLa_10bin, DNAse)
save(HeLa_tiled, file = "~/Human_Datasets/HeLa/HeLa_after_DNAse.Rda")

load("~/Human_Datasets/HeLa/HeLa_after_DNAse.Rda")
# Adding STARR-seq annotations
STARR_t1 <- read.csv("~/Human_Datasets/HeLa/STARR_seq/peaks_inhibitor_correctedEnrichment4_supp.table3.tsv",
  sep = "\t")

STARR_t2 <- read.csv("/home/jw18713/Human_Datasets/HeLa/STARR_seq/shortlisted_regions_inhibitor_supp.table3.tsv",
  sep = "\t")

STARR_1 <- GRanges(seqnames = STARR_t1$seqnames, ranges = IRanges(
  start = STARR_t1$start,
  end = STARR_t1$end))

STARR_2 <- GRanges(seqnames = STARR_t2$seqnames, ranges = IRanges(
  start = STARR_t2$start,
  end = STARR_t2$end))

STARR <- c(STARR_1, STARR_2)

# STARR-seq liftOver
STARR <- unlist(liftOver(STARR, chain))

# Assigning a 1 to STARR-seq overlaps
overlaps <- findOverlaps(HeLa_tiled,STARR)
scores <- rep(0, length(HeLa_tiled))
scores[as.vector(queryHits(overlaps))] <- 1
HeLa_tiled$STARR_seq_binary <- scores

# Saving the lifted over STARR-seq dataset
save(HeLa_tiled, file = "~/Human_Datasets/HeLa/DMRcaller_with_STARR_DNAse_lenient.Rda")

# Normalising the columns by min of max normalisation
for (c in seq_along(mcols(HeLa_tiled))){
  max_value = max(mcols(HeLa_tiled)[c][,1])
  min_value = min(mcols(HeLa_tiled)[c][,1])
  mcols(HeLa_tiled)[c][,1] <- (mcols(HeLa_tiled)[c][,1] - min_value) / (max_value - min_value)
}

save(HeLa_tiled, file = "~/Human_Datasets/HeLa/DMRcaller_norm_tiled_DNAse_lenient.Rda")
