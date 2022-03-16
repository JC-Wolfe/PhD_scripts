
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# HeLa matrix
HeLa_cmat <- matrix(0,2,3)
colnames(HeLa_cmat) <- c("Accuracy", "Precision", "Recall")
rownames(HeLa_cmat) <- c("XAI", "NN")

HeLa_cmat[1,1] <- 0.42188
HeLa_cmat[1,2] <- 0.51094
HeLa_cmat[1,3] <- 0.58333

HeLa_cmat[2,1] <- 0.52503
HeLa_cmat[2,2] <- 0.51126
HeLa_cmat[2,3] <- 0.58897

# Plotting HeLa
pdf(file = "~/Human_Datasets/Plots/Figure1/Extra_Performance/HeLa_lenient_lenient_barplots.pdf",
width = 6, height = 6, pointsize = 14)
par(mar=c(4,3,4,5.5), cex = 1.2)
par(xpd = NA)
    barplot(HeLa_cmat, col = cbbPalette[2:3], beside = T, ylim = c(0,1),
            main = "HeLa Lenient:Lenient Statistics"
            )
    legend(10,0.75,rownames(HeLa_cmat), fill=cbbPalette[2:3])
dev.off()
