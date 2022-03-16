
library(GenomicRanges)
library(rtracklayer)

setwd("~/P1_Full_Run/S2_Model_outputs/BG3/XAI")

print("BG3 XAI Start")

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

BG3_XAI <- GRanges()

for (i in dir()){
  BG3_XAI <- c(BG3_XAI, read_LG_output(i))
}

print("BG3 XAI End")

save(BG3_XAI, file = "~/P1_Full_Run/S2_Model_outputs/BG3_XAI.Rda")
rm(list=ls())
