install.packages("PRROC")
library(PRROC)
# Load the required packages
require(ggplot2)
require(ggseqlogo)

load("~/Human_Datasets/H9/Processed_datasets/H9_XAI.Rda")

fg1 <- H9_XAI$Conf_Perc_1[H9_XAI$STARR_seq_binary == 1]
bg1 <- H9_XAI$Conf_Perc_1[H9_XAI$STARR_seq_binary == 0]

# H9 ROC Curve
H9_roc <- roc.curve(scores.class0 = fg1, scores.class1 = bg1, curve = T)

# H9 PR Curve
H9_pr <- pr.curve(scores.class0 = fg1, scores.class1 = bg1, curve = T)

png("~/Human_Datasets/H9/Plots/AUPRC/ROC.png", width=500, height=500)
plot(H9_roc, main = "H9 ROC")
dev.off()
png("~/Human_Datasets/H9/Plots/AUPRC/PR.png", width=500, height=500)
plot(H9_pr, main = "H9 PRC")
dev.off()
