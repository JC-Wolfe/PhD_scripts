library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# S2 matrix
S2_cmat <- matrix(0,2,3)
colnames(S2_cmat) <- c("Accuracy", "Precision", "Recall")
rownames(S2_cmat) <- c("XAI", "NN")

S2_cmat[1,1] <- 0.8260717
S2_cmat[1,2] <- 0.5739867
S2_cmat[1,3] <- 0.7134783

S2_cmat[2,1] <- 0.79425
S2_cmat[2,2] <- 0.5715117
S2_cmat[2,3] <- 0.7397067

# BG3 Matrix

BG3_cmat <- matrix(0,2,3)
colnames(BG3_cmat) <- c("Accuracy", "Precision", "Recall")
rownames(BG3_cmat) <- c("XAI", "NN")

BG3_cmat[1,1] <- 0.36621
BG3_cmat[1,2] <- 0.5064617
BG3_cmat[1,3] <- 0.6169067

BG3_cmat[2,1] <- 0.2670167
BG3_cmat[2,2] <- 0.5073833
BG3_cmat[2,3] <- 0.612255

# Plotting BG3
pdf(file = "~/Review_Paper_Plots/arch_prot_rules/plots/Barplots/S2_Barplots(LG).pdf",
width = 12, height = 6, pointsize = 14)
par(mar=c(4,3,4,5.5), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)
    barplot(S2_cmat, col = cbbPalette[2:3], beside = T, ylim = c(0,1),
            main = "S2 Confusion Matrix Statistics"
            )
    barplot(BG3_cmat, col = cbbPalette[2:3], beside = T, ylim = c(0,1),
            main = "BG3 Confusion Matrix Statistics"
            )
    legend(10,0.75,rownames(S2_cmat), fill=cbbPalette[2:3])
dev.off()
