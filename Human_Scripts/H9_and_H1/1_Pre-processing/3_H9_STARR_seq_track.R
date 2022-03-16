library(GenomicRanges)
library(rtracklayer)
library(genomation)

load("~/Human_Datasets/H9/DMRcaller_norm_tiled.Rda")

H9_STARR_seq <- H9_tiled[H9_tiled$STARR_seq_binary == 1]
H9_STARR_seq <- reduce(H9_STARR_seq)

save(H9_STARR_seq, file = "~/Human_Datasets/H9/H9_STARR_seq.Rda")
