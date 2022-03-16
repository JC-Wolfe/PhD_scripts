install.packages("PRROC")
library(PRROC)
# Load the required packages
require(ggplot2)
require(ggseqlogo)

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


fg1 <- DvHK_XAI$Conf_Perc_1[DvHK_XAI$Developmental == 1]
bg1 <- DvHK_XAI$Conf_Perc_1[DvHK_XAI$Developmental == 0]

# DvHK ROC Curve
DvHK_roc <- roc.curve(scores.class0 = fg1, scores.class1 = bg1, curve = T)

# DvHK PR Curve
DvHK_pr <- pr.curve(scores.class0 = fg1, scores.class1 = bg1, curve = T)



png("/home/jw18713/Project1/Paper_Plots/FigureS1/DvHk_roc.png", width=500, height=500)
plot(DvHK_roc, main = "Developmental vs Housekeeping ROC")
dev.off()
png("/home/jw18713/Project1/Paper_Plots/FigureS1/DvHk_pr.png", width=500, height=500)
plot(DvHK_pr, main = "Developmental vs Housekeeping PRC")
dev.off()
