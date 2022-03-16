
library(GenomicRanges)
library(rtracklayer)

setwd("~/Human_Datasets/H9/csv_files/XAI_Processed")

print("H9 XAI Start")

read_LG_output <- function(x){
  out_df <- read.delim(x, sep = "\t", stringsAsFactors = F)
  out_gr <- GRanges(seqnames = out_df$seqnames,
    ranges = IRanges(start = out_df$start, end = out_df$end),
    strand = out_df$strand)
  out_gr$STARR_seq_binary <- out_df$STARR_seq_binary
  out_gr$Predicted_Class <- out_df$Predicted.Class
  out_gr$Conf_Perc_1 <- out_df$Conf..Perc.1
  out_gr$Conf_Perc_0 <- out_df$Conf..Perc.0
  out_gr$Errors <- out_df$Errors
  return(out_gr)
}

H9_XAI <- GRanges()

for (i in dir()){
  H9_XAI <- c(H9_XAI, read_LG_output(i))
}

print("H9 XAI End")

save(H9_XAI, file = "~/Human_Datasets/H9/Processed_datasets/H9_XAI.Rda")
rm(list=ls())
