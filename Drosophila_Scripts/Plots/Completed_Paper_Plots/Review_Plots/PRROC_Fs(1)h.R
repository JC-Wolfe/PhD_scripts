install.packages("PRROC")
library(PRROC)
# Load the required packages
require(ggplot2)
require(ggseqlogo)

Arch_bulk <- read.csv("~/Review_Paper_Plots/arch_prot_rules/data/No_fs(1)h_bulk.tsv", sep = "\t")

fg1 <- Arch_bulk$Conf..Perc.1[Arch_bulk$IsCommon == 1]
bg1 <- Arch_bulk$Conf..Perc.1[Arch_bulk$IsCommon == 0]

# Arch ROC Curve
Arch_roc <- roc.curve(scores.class0 = fg1, scores.class1 = bg1, curve = T)

# Arch PR Curve
Arch_pr <- pr.curve(scores.class0 = fg1, scores.class1 = bg1, curve = T)


png("~/Review_Paper_Plots/arch_prot_rules/plots/PRROC/Fs(1)h_ROC.png", width=500, height=500)
plot(Arch_roc, main = "Chromatin Features ROC")
dev.off()
png("~/Review_Paper_Plots/arch_prot_rules/plots/PRROC/Fs(1)h_PRC.png", width=500, height=500)
plot(Arch_pr, main = "Chromatin Features PRC")
dev.off()
