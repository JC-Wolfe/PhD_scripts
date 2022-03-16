library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)

load("~/Human_Datasets/H1/Processed_datasets/H1_full(widths).Rda")
load("~/Human_Datasets/H9/Processed_datasets/H9_full(widths).Rda")

H9_putative <- reduce(H9_full[H9_full$classification == "Putative"])
H9_common <- reduce(H9_full[H9_full$classification == "Common"])

H1_common <- reduce(H1_full[H1_full$classification == "Common"])
H1_specific <- reduce(H1_full[H1_full$classification == "H1 Specific"])

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# H9_putative
pdf("~/Human_Datasets/Plots/Figure2/Figure2C_H9_Putative(width).pdf",
width = 9, height = 8, pointsize = 14)
par(mar=c(4,6.5,4,4), cex = 1.2)
par(xpd = NA)

hist(log10(width(H9_putative)), axes=F, main="H9 Putative Enhancer Widths",
  xlab="Width (bp)", col=cbbPalette[3], cex.main=2, cex.lab=1.5, xlim = c(2,5),
  ylim = c(0,100000), breaks = 12)
axis(1, at=c(log10(100), log10(1000), log10(10000), log10(100000)),
    labels=c("100", "1000", "10000", "100000"))
axis(2, at=seq(0,100000, by = 25000),
  labels = c("0", "25k", "50k", "75k", "100k"), las = 1)

dev.off()

# H9_common
pdf("~/Human_Datasets/Plots/Figure2/Figure2C_H9_Common(width).pdf",
width = 9, height = 8, pointsize = 14)
par(mar=c(4,6.5,4,4), cex = 1.2)
par(xpd = NA)

hist(log10(width(H9_common)), axes=F, main="H9 Common Enhancer Widths",
  xlab="Width (bp)", col=cbbPalette[3], cex.main=2, cex.lab=1.5, xlim = c(2,5),
  ylim = c(0,12000), breaks = 12)
axis(1, at=c(log10(100), log10(1000), log10(10000), log10(100000)),
    labels=c("100", "1000", "10000", "100000"))
axis(2, at=seq(0,12000, by = 3000),
  labels = c("0", "3000", "6000", "9000", "12000"), las = 1)

dev.off()

# H1_common
pdf("~/Human_Datasets/Plots/Figure2/Figure2C_H1_Common(width).pdf",
width = 9, height = 8, pointsize = 14)
par(mar=c(4,6.5,4,4), cex = 1.2)
par(xpd = NA)

hist(log10(width(H1_common)), axes=F, main="H1 Common Enhancer Widths",
  xlab="Width (bp)", col=cbbPalette[3], cex.main=2, cex.lab=1.5, xlim = c(2,5),
  ylim = c(0,75000), breaks = 12)
axis(1, at=c(log10(100), log10(1000), log10(10000), log10(100000)),
    labels=c("100", "1000", "10000", "100000"))
axis(2, at=seq(0,75000, by = 25000),
  labels = c("0", "25k", "50k", "75k"), las = 1)

dev.off()

# H1_specific
pdf("~/Human_Datasets/Plots/Figure2/Figure2C_H1_Specific(width).pdf",
width = 9, height = 8, pointsize = 14)
par(mar=c(4,6.5,4,4), cex = 1.2)
par(xpd = NA)

hist(log10(width(H1_specific)), axes=F, main="H1 Specific Enhancer Widths",
  xlab="Width (bp)", col=cbbPalette[3], cex.main=2, cex.lab=1.5, xlim = c(2,5),
  ylim = c(0,50000), breaks = 12)
axis(1, at=c(log10(100), log10(1000), log10(10000), log10(100000)),
    labels=c("100", "1000", "10000", "100000"))
axis(2, at=seq(0,50000, by = 10000),
  labels = c("0", "10k", "20k", "30k", "40k", "50k"), las = 1)

dev.off()
