library(GenomicRanges)
library(rtracklayer)

# Putative
load("~/Paper_replots/Objects/Class_split_glists/BG3_putative_list.Rda")
load("~/Paper_replots/Objects/Class_split_glists/S2_putative_list.Rda")

# Common
load("~/Paper_replots/Objects/Class_split_glists/BG3_shared_list.Rda")
load("~/Paper_replots/Objects/Class_split_glists/S2_shared_list.Rda")

# STARR-seq
load("~/Paper_replots/Objects/Class_split_glists/BG3_starr_res_list.Rda")
load("~/Paper_replots/Objects/Class_split_glists/S2_starr_res_list.Rda")

# Background
load("~/Paper_replots/Objects/Class_split_glists/BG3_bg_list.Rda")
load("~/Paper_replots/Objects/Class_split_glists/S2_bg_list.Rda")

Neither_matrix <- matrix(0, 2, 4, dimnames = list(c("BG3", "S2"),
  c("Putative", "Common", "STARR-seq only", "Whole Genome")))

Neither_matrix[1,1] <- length(BG3_putative[[4]])
Neither_matrix[1,2] <- length(BG3_shared[[4]])
Neither_matrix[1,3] <- length(BG3_starr_res[[4]])
Neither_matrix[1,4] <- length(BG3_bg[[4]])

Neither_matrix[2,1] <- length(S2_putative[[4]])
Neither_matrix[2,2] <- length(S2_shared[[4]])
Neither_matrix[2,3] <- length(S2_starr_res[[4]])
Neither_matrix[2,4] <- length(S2_bg[[4]])

save(Neither_matrix, file = "~/Paper_replots/Objects/Neither_matrix.Rda")
