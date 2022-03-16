library(pROC)

read_LG_output <- function(x){
  out_df <- read.delim(x, sep = "\t", stringsAsFactors = F)
  out_gr <- GRanges(seqnames = out_df$seqnames,
    ranges = IRanges(start = out_df$start, end = out_df$end),
    strand = out_df$strand)
  out_gr$Developmental <- out_df$Developmental
  out_gr$Predicted_Class <- out_df$Predicted.Class
  out_gr$Conf_Perc_1 <- out_df$Conf..Perc.1
  out_gr$Conf_Perc_0 <- out_df$Conf..Perc.0
  out_gr$Errors <- out_df$Errors
  return(out_gr)
}

DvHK_XAI <- read_LG_output("dev_hk_results.tsv")

v1 <- roc(DvHK_XAI$Developmental, DvHK_XAI$Conf_Perc_1)

pdf("/home/jw18713/Project1/Paper_Plots/FigureS1/DvHk_B&W_roc.pdf")
plot(v1, main = "Developmental vs Housekeeping ROC")
dev.off()
