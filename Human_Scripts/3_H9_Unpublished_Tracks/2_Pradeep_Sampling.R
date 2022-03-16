
library(GenomicRanges)
library(rtracklayer)
library(genomation)
options(scipen=999)

load("~/Human_Datasets/H9/H9_Pradeep_NAs.Rda")

H9_tiled <- H9_tiled[!is.na(H9_tiled$STARR_seq_binary)]

chr_ratios <- width(seqnames(H9_tiled))/sum(width(seqnames(H9_tiled)))

positive_bins <- H9_tiled[H9_tiled$STARR_seq_binary == 1]
pos_num <- round(chr_ratios * length(positive_bins)/length(H9_tiled) * 500000)
negative_bins <- H9_tiled[H9_tiled$STARR_seq_binary == 0]
neg_num <- round(chr_ratios * length(negative_bins)/length(H9_tiled) * 500000)

pos_split <- split(positive_bins, seqnames(positive_bins))
neg_split <- split(negative_bins, seqnames(negative_bins))

# set seed so sampled regions are repeatable
set.seed(1809)

# So lets get those negative samples out then
H9_input_500k <- GRanges()
H9_unseen_500k <- GRanges()

for (i in seq_along(neg_split)){
  sbins <- sample(neg_split[[i]], size = neg_num[i], replace = F)
  remaining <- neg_split[[i]][neg_split[[i]] %over% setdiff(neg_split[[i]], sbins)]
  remaining_sample <- sample(remaining, size = neg_num[i], replace = F)
  H9_input_500k <- c(H9_input_500k, sbins)
  H9_unseen_500k <- c(H9_unseen_500k, remaining_sample)
}

# And the positive
for (i in seq_along(pos_split)){
  sbins <- sample(pos_split[[i]], size = pos_num[i], replace = F)
  remaining <- pos_split[[i]][pos_split[[i]] %over% setdiff(pos_split[[i]], sbins)]
  remaining_sample <- sample(remaining, size = pos_num[i], replace = F)
  H9_input_500k <- c(H9_input_500k, sbins)
  H9_unseen_500k <- c(H9_unseen_500k, remaining_sample)
}

# And export the files we need
save(H9_input_500k, file = "~/Human_Datasets/H9/Pradeep/H9_inc_Pradeep.Rda")
write.csv(H9_input_500k, file = "~/Human_Datasets/H9/Pradeep/H9_inc_Pradeep.csv",
  row.names = F)

# And the holdout files
save(H9_unseen_500k, file = "~/Human_Datasets/H9/Pradeep/H9_inc_Pradeep_holdout.Rda")
write.csv(H9_unseen_500k, file = "~/Human_Datasets/H9/Pradeep/H9_inc_Pradeep_holdout.csv",
  row.names = F)
