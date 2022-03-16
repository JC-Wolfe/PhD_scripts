library(GenomicRanges)
library(rtracklayer)
library(gplots)
library(lattice)

setwd("~/Pre_viva/Hi-C")
load("contact_Rda/H9_putative.Rda")
load("contact_Rda/H9_common.Rda")
load("contact_Rda/H9_starr_only.Rda")
load("contact_Rda/H9_neither.Rda")
load("contact_Rda/H9_background.Rda")

H9_mat <- matrix(0,4,4)
rownames(H9_mat) <- c("Putative", "Common", "STARR-seq", "Background")
colnames(H9_mat) <- c("Distal & Proximal", "Distal Only", "Proximal Only", "Neither")

H9_mat[1,1] <- length(H9_common[[1]])
H9_mat[1,2] <- length(H9_common[[2]])
H9_mat[1,3] <- length(H9_common[[3]])
H9_mat[1,4] <- length(H9_common[[4]])

H9_mat[2,1] <- length(H9_putative[[1]])
H9_mat[2,2] <- length(H9_putative[[2]])
H9_mat[2,3] <- length(H9_putative[[3]])
H9_mat[2,4] <- length(H9_putative[[4]])

H9_mat[3,1] <- length(H9_starr_only[[1]])
H9_mat[3,2] <- length(H9_starr_only[[2]])
H9_mat[3,3] <- length(H9_starr_only[[3]])
H9_mat[3,4] <- length(H9_starr_only[[4]])

H9_mat[4,1] <- length(H9_background[[1]])
H9_mat[4,2] <- length(H9_background[[2]])
H9_mat[4,3] <- length(H9_background[[3]])
H9_mat[4,4] <- length(H9_background[[4]])

H9_bar <- t(H9_mat)
H9_bar_aves <- H9_bar

H9_bar_aves[,1] <- H9_bar[,1]/sum(H9_bar[,1])
H9_bar_aves[,2] <- H9_bar[,2]/sum(H9_bar[,2])
H9_bar_aves[,3] <- H9_bar[,3]/sum(H9_bar[,3])
H9_bar_aves[,4] <- H9_bar[,4]/sum(H9_bar[,4])

cs <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

setwd("~/Pre_viva/Hi-C/Plots")

pdf(paste0("Figure3A1.pdf"),
width = 11, height = 8, pointsize = 14)
par(mar=c(4,4,4,10.5), cex = 1.2)
par(xpd = NA)

barplot(H9_bar_aves, col=cs, xlab="Classification",
  ylab="Percentage of Results With Contacts",
  main="H9 Contact Comparisons", yaxt = "none")
axis(2, seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10), "%"), las=2)
legend(5.25,0.75,rownames(H9_bar_aves), fill=cs, bty = "n")

dev.off()

H9_hmap <- H9_bar
H9_hmap_aves <- H9_bar_aves


H9_hmap_comp <- H9_hmap_aves[,1:3]
H9_hmap_comp[1,] <- log2(H9_hmap_comp[1,]/H9_hmap_aves[1,4])
H9_hmap_comp[2,] <- log2(H9_hmap_comp[2,]/H9_hmap_aves[2,4])
H9_hmap_comp[3,] <- log2(H9_hmap_comp[3,]/H9_hmap_aves[3,4])
H9_hmap_comp[4,] <- log2(H9_hmap_comp[4,]/H9_hmap_aves[4,4])

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
custom_at <- seq(-1.2, 1.2, by = (2.4/30))
cols_contrast <- rev(c(colorRampPalette(c(cbbPalette[7], "white", cbbPalette[6]))(30)))
library(gridExtra)

v1 <- levelplot(H9_hmap_comp,
  at = custom_at,
  main = list(label = "Predicted Region Contact Comparisons H9", cex = 2.4),
  xlab = list(label = "log2 Enrichment Difference", cex = 1.8),
  ylab = list(label = "Predicted By", cex = 1.8),
  scales=list(x=list(rot=90, cex=1.2), y=list(cex=1.2)),
  col.regions = cols_contrast)

pdf(paste0("Figure3A2.pdf"),
width = 22, height = 8, pointsize = 14)
par(mar=c(4,4,4,10.5), cex = 1.2)
grid.arrange(v1, nrow = 1)
dev.off()
