
library(genomation)
library(GenomicRanges)
library(rtracklayer)

setwd("~/Human_Datasets/H9/Extra_tracks")
load("Processed_classifications/putative(widths).Rda")
load("Processed_classifications/common(widths).Rda")
load("Processed_classifications/background(widths).Rda")

pmat <- as.matrix(mcols(putative))
bgmat <- as.matrix(mcols(ET_BG))
cmat <- as.matrix(mcols(common))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

setwd("~/Human_Datasets/Plots/Figure6_boxplots/toproom")

# PolII
pdf("PolII_boxplot(widths).pdf",
  width = 8, height = 6, pointsize = 14)
par(mar = c(4,5.5,4,8), xpd = T)

boxplot(pmat[,1], cmat[,1], bgmat[,1],
  names = c("Putative", "Common", "Background"),
  ylim = c(0, 1.4),
  outline = F, col = cbbPalette[2:4], main = colnames(bgmat)[1])
dev.off()

# PRO_seq
pdf("PRO_seq_boxplot(widths).pdf",
  width = 8, height = 6, pointsize = 14)
par(mar = c(4,5.5,4,8), xpd = T)

boxplot(pmat[,2], cmat[,2], bgmat[,2],
  names = c("Putative", "Common", "Background"),
  ylim = c(0, 3.5),
  outline = F, col = cbbPalette[2:4], main = colnames(bgmat)[2])
dev.off()

# Rad21
pdf("Rad21_boxplot(widths).pdf",
  width = 8, height = 6, pointsize = 14)
par(mar = c(4,5.5,4,8), xpd = T)

boxplot(pmat[,3], cmat[,3], bgmat[,3],
  names = c("Putative", "Common", "Background"),
  ylim = c(0, 12),
  outline = F, col = cbbPalette[2:4], main = colnames(bgmat)[3])
dev.off()

# MED1
pdf("MED1_boxplot(widths).pdf",
  width = 8, height = 6, pointsize = 14)
par(mar = c(4,5.5,4,8), xpd = T)

boxplot(pmat[,4], cmat[,4], bgmat[,4],
  names = c("Putative", "Common", "Background"),
  ylim = c(0, 1.4),
  outline = F, col = cbbPalette[2:4], main = colnames(bgmat)[4])
dev.off()
