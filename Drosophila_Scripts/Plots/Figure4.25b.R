library(pROC)

Arch <- read.csv("~/Review_Paper_Plots/arch_prot_rules/data/No_fs(1)h_bulk.tsv", sep = "\t")

v1 <- roc(Arch$IsCommon, Arch$Conf..Perc.1)

pdf("~/Review_Paper_Plots/arch_prot_rules/plots/PRROC/No_Fs(1)h_B&W_roc.pdf")
plot(v1, main = "Common vs Putative ROC")
dev.off()
