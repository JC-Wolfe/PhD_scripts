
library(GenomicRanges)
library(rtracklayer)
library(genomation)
options(scipen=999)

load("~/Human_Datasets/H1/DMRcaller_norm_tiled.Rda")

# I need to figure out the ratios of chromosomes to be included in the training and testing data
# I'm specifying 1:6 because there is no chrY data for S2 cells

chr_ratios <- width(seqnames(H1_tiled))/sum(width(seqnames(H1_tiled)))

# Now that I have the ratios I need to figure out how many samples I need from each chromosome to
# proportionally represent them in the 2 million data points I will be using. I'm rounding to be
# multiples of 100
positive_bins <- H1_tiled[H1_tiled$STARR_seq_binary == 1]
pos_num <- round(chr_ratios * length(positive_bins)/length(H1_tiled) * 500000)
negative_bins <- H1_tiled[H1_tiled$STARR_seq_binary == 0]
neg_num <- round(chr_ratios * length(negative_bins)/length(H1_tiled) * 500000)

pos_split <- split(positive_bins, seqnames(positive_bins))
neg_split <- split(negative_bins, seqnames(negative_bins))

# set seed so sampled regions are repeatable
set.seed(1509)

# So lets get those negative samples out then
H1_input_500k <- GRanges()
H1_unseen_500k <- GRanges()

for (i in seq_along(neg_split)){
  sbins <- sample(neg_split[[i]], size = neg_num[i], replace = F)
  remaining <- neg_split[[i]][neg_split[[i]] %over% setdiff(neg_split[[i]], sbins)]
  remaining_sample <- sample(remaining, size = neg_num[i], replace = F)
  H1_input_500k <- c(H1_input_500k, sbins)
  H1_unseen_500k <- c(H1_unseen_500k, remaining_sample)
}

# And the positive
for (i in seq_along(pos_split)){
  sbins <- sample(pos_split[[i]], size = pos_num[i], replace = F)
  remaining <- pos_split[[i]][pos_split[[i]] %over% setdiff(pos_split[[i]], sbins)]
  remaining_sample <- sample(remaining, size = pos_num[i], replace = F)
  H1_input_500k <- c(H1_input_500k, sbins)
  H1_unseen_500k <- c(H1_unseen_500k, remaining_sample)
}


# And export the files we need
save(H1_input_500k, file = "~/Human_Datasets/H1/DMR_method_input_500k.Rda")
write.csv(H1_input_500k, file = "~/Human_Datasets/H1/DMR_method_input_500k.csv",
  row.names = F)

# And the holdout files
save(H1_unseen_500k, file = "~/Human_Datasets/H1/DMR_method_holdout_500k.Rda")
write.csv(H1_unseen_500k, file = "~/Human_Datasets/H1/DMR_method_holdout_500k.csv",
  row.names = F)


# Now to check the method with NAs
H1_input_NAs <- H1_input_500k
H1_unseen_NAs <- H1_unseen_500k

zero_to_NA <- function(x){
  x[x == 0] <- NA
  return(x)
}

NA_cols <- apply(mcols(H1_input_NAs)[1:32], 2, zero_to_NA)
mcols(H1_input_NAs)[1:32] <- NA_cols

NA_unseen <- apply(mcols(H1_unseen_NAs)[1:32], 2, zero_to_NA)
mcols(H1_unseen_NAs)[1:32] <- NA_unseen


# And export the NA files
save(H1_input_NAs, file = "~/Human_Datasets/H1/DMR_method_input_NAs.Rda")
write.csv(H1_input_NAs, file = "~/Human_Datasets/H1/DMR_method_input_NAs.csv",
  row.names = F)

# And the holdout NA files
save(H1_unseen_NAs, file = "~/Human_Datasets/H1/DMR_method_holdout_NAs.Rda")
write.csv(H1_unseen_NAs, file = "~/Human_Datasets/H1/DMR_method_holdout_NAs.csv",
  row.names = F)
