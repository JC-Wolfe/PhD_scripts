library(genomation)
library(GenomicRanges)
library(rtracklayer)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# H9 matrix
H9_cmat <- matrix(0,2,3)
colnames(H9_cmat) <- c("Accuracy", "Precision", "Recall")
rownames(H9_cmat) <- c("XAI", "NN")

H9_cmat[1,1] <- 0.72535
H9_cmat[1,2] <- 0.50869
H9_cmat[1,3] <- 0.73059

H9_cmat[2,1] <- 0.78653
H9_cmat[2,2] <- 0.51071
H9_cmat[2,3] <- 0.74007

# Plotting H9
pdf(file = "~/Human_Datasets/Plots/Figure1/H9_confusion_matrix_barplots.pdf",
width = 6, height = 6, pointsize = 14)
par(mar=c(4,3,4,5.5), cex = 1.2)
par(xpd = NA)
    barplot(H9_cmat, col = cbbPalette[2:3], beside = T, ylim = c(0,1),
            main = "H9 Confusion Matrix Statistics"
            )
    legend(10,0.75,rownames(H9_cmat), fill=cbbPalette[2:3])
dev.off()
