
library(GenomicRanges)
library(rtracklayer)
library(genomation)
options(scipen=999)

load("~/Human_Datasets/HeLa/DMRcaller_norm_tiled_DNAse_strict.Rda")

zero_to_NA <- function(x){
  x[x == 0] <- NA
  return(x)
}

NA_cols <- apply(mcols(HeLa_tiled)[1:12], 2, zero_to_NA)
mcols(HeLa_tiled)[1:12] <- NA_cols

chr_ratios <- width(seqnames(HeLa_tiled))/sum(width(seqnames(HeLa_tiled)))

positive_bins <- HeLa_tiled[HeLa_tiled$STARR_seq_binary == 1]
pos_num <- round(chr_ratios * length(positive_bins)/length(HeLa_tiled) * 500000)
negative_bins <- HeLa_tiled[HeLa_tiled$STARR_seq_binary == 0]
neg_num <- round(chr_ratios * length(negative_bins)/length(HeLa_tiled) * 500000)

pos_split <- split(positive_bins, seqnames(positive_bins))
neg_split <- split(negative_bins, seqnames(negative_bins))

# set seed so sampled regions are repeatable
set.seed(1809)

# So lets get those negative samples out then
HeLa_input_500k <- GRanges()
HeLa_unseen_500k <- GRanges()

for (i in seq_along(neg_split)){
  sbins <- sample(neg_split[[i]], size = neg_num[i], replace = F)
  remaining <- neg_split[[i]][neg_split[[i]] %over% setdiff(neg_split[[i]], sbins)]
  remaining_sample <- sample(remaining, size = neg_num[i], replace = F)
  HeLa_input_500k <- c(HeLa_input_500k, sbins)
  HeLa_unseen_500k <- c(HeLa_unseen_500k, remaining_sample)
}

# And the positive
for (i in seq_along(pos_split)){
  sbins <- sample(pos_split[[i]], size = pos_num[i], replace = F)
  remaining <- pos_split[[i]][pos_split[[i]] %over% setdiff(pos_split[[i]], sbins)]
  remaining_sample <- sample(remaining, size = pos_num[i], replace = F)
  HeLa_input_500k <- c(HeLa_input_500k, sbins)
  HeLa_unseen_500k <- c(HeLa_unseen_500k, remaining_sample)
}

# And export the files we need
save(HeLa_input_500k, file = "~/Human_Datasets/HeLa/HeLa_NAs_input_strict.Rda")
write.csv(HeLa_input_500k, file = "~/Human_Datasets/HeLa/HeLa_NAs_input_strict.csv",
  row.names = F)

# And the holdout files
save(HeLa_unseen_500k, file = "~/Human_Datasets/HeLa/HeLa_NAs_holdout_strict.Rda")
write.csv(HeLa_unseen_500k, file = "~/Human_Datasets/HeLa/HeLa_NAs_holdout_strict.csv",
  row.names = F)
