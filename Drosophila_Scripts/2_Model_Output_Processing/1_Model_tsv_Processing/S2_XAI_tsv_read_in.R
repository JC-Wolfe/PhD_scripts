
library(GenomicRanges)
library(rtracklayer)

setwd("~/P1_Full_Run/Model_outputs/S2_tsv/XAI")

print("S2 XAI Start")

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

S2_XAI <- GRanges()

for (i in dir()){
  S2_XAI <- c(S2_XAI, read_LG_output(i))
}

print("S2 XAI End")

save(S2_XAI, file = "~/P1_Full_Run/Model_outputs/S2_XAI.Rda")
rm(list=ls())
