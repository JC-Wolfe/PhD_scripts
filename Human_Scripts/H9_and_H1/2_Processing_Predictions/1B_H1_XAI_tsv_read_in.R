
library(GenomicRanges)
library(rtracklayer)

setwd("~/Human_Datasets/H1/csv_files/XAI_processed")

print("H1 XAI Start")

read_LG_output <- function(x){
  out_df <- read.delim(x, sep = "\t", stringsAsFactors = F)
  out_gr <- GRanges(seqnames = out_df$seqnames,
    ranges = IRanges(start = out_df$start, end = out_df$end),
    strand = out_df$strand)
  # There is no STARR-seq in H1 cells
  # out_gr$STARR_seq_binary <- out_df$STARR_seq_binary
  out_gr$Predicted_Class <- out_df$Predicted.Class
  out_gr$Conf_Perc_1 <- out_df$Conf..Perc.1
  out_gr$Conf_Perc_0 <- out_df$Conf..Perc.0
  out_gr$Errors <- out_df$Errors
  return(out_gr)
}

H1_XAI <- GRanges()

for (i in dir()){
  H1_XAI <- c(H1_XAI, read_LG_output(i))
}

print("H1 XAI End")

save(H1_XAI, file = "~/Human_Datasets/H1/Processed_datasets/H1_XAI.Rda")

H1_predicted <- H1_XAI[H1_XAI$Predicted_Class == 1]

save(H1_predicted, file = "~/Human_Datasets/H1/Processed_datasets/H1_predicted.Rda")
