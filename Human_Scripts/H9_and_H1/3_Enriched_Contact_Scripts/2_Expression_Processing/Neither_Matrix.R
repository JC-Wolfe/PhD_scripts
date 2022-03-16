library(GenomicRanges)
library(rtracklayer)
library(gplots)
library(lattice)

setwd("~/Pre_viva/Hi-C")
load("contact_Rda/H9_putative.Rda")
load("contact_Rda/H9_common.Rda")
load("contact_Rda/H9_starr_only.Rda")
load("contact_Rda/H9_background.Rda")

Neither_matrix <- matrix(0, 1, 4, dimnames = list(c("H9"),
  c("Putative", "Common", "STARR-seq Only", "Whole Genome")))

Neither_matrix[1,1] <- length(H9_putative[[4]])
Neither_matrix[1,2] <- length(H9_common[[4]])
Neither_matrix[1,3] <- length(H9_starr_only[[4]])
Neither_matrix[1,4] <- length(H9_background[[4]])

setwd("~/Pre_viva/Hi-C/grange_lists")
save(Neither_matrix, file = "Neither_matrix.Rda")
