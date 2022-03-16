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
Pradeep_tiled <- empty_tiled_genome(Hsapiens, binsize = 100, chr_range = seq(1,23))

# Next I need an empty tiled genome set at every 10 base pairs to average the
# reads over
Pradeep_10bin <- empty_tiled_genome(Hsapiens, binsize = 10, chr_range = seq(1,23))

# Now I need to populate the mcols with track information
setwd("~/Human_Datasets/Pradeep/Pradeep")

files <- dir()

for (i in seq_along(files)){
    t_read <- read.csv(files[i], sep = "\t", header = F)
    track <- GRanges(seqnames = t_read[,1], ranges = IRanges(start = t_read[,2], end = t_read[,3]))
    track$score <- t_read[,4]
    mcols(Pradeep_tiled)[length(mcols(Pradeep_tiled)) + 1] <- track_calculator(Pradeep_10bin,
      track)
    sname <- unlist(strsplit(files[i], "_"))
    sname <- sname[!sname == "treat" & !sname == "pileup.bdg"]
    fname <- paste(sname, collapse = "_")
    names(mcols(Pradeep_tiled))[length(mcols(Pradeep_tiled))] <- fname
}

# Normalising the columns by min of max normalisation
for (c in seq_along(mcols(Pradeep_tiled))){
  max_value = max(mcols(Pradeep_tiled)[c][,1])
  min_value = min(mcols(Pradeep_tiled)[c][,1])
  mcols(Pradeep_tiled)[c][,1] <- (mcols(Pradeep_tiled)[c][,1] - min_value) / (max_value - min_value)
}

# Writing NAs where appropriate
zero_to_NA <- function(x){
  x[x == 0] <- NA
  return(x)
}

NA_cols <- apply(mcols(Pradeep_tiled)[1:10], 2, zero_to_NA)
mcols(Pradeep_tiled)[1:10] <- NA_cols

# Saving the histone modifications and clearing memory
save(Pradeep_tiled, file = "~/Human_Datasets/H9/Pradeep/Pradeep_Tracks.Rda")

# Loading in H9
load("~/Human_Datasets/H9/DMRcaller_norm_tiled_DNAse.Rda")

# Setting 0 to NA in H9
NA_cols <- apply(mcols(H9_tiled)[1:32], 2, zero_to_NA)
mcols(H9_tiled)[1:32] <- NA_cols

H9_tiled$H3K122_ab_old <- Pradeep_tiled$H3K122_ab_old

H9_tiled$H3K122_abc_new <- Pradeep_tiled$H3K122_abc_new

H9_tiled$H3K122ac_combined <- Pradeep_tiled$H3K122ac_combined

H9_tiled$H3K4me2_1000 <- Pradeep_tiled$H3K4me2_1000

H9_tiled$H3K4me2_250 <- Pradeep_tiled$H3K4me2_250

H9_tiled$H3K4me2_500 <- Pradeep_tiled$H3K4me2_500

H9_tiled$H3K4me2_combined <- Pradeep_tiled$H3K4me2_combined

H9_tiled$H4K16ac_aug <- Pradeep_tiled$H4K16ac_aug

H9_tiled$H4K16ac_combined <- Pradeep_tiled$H4K16ac_combined

H9_tiled$H4K16ac_sep <- Pradeep_tiled$H4K16ac_sep

save(H9_tiled, file = "~/Human_Datasets/H9/H9_Pradeep_NAs.Rda")
