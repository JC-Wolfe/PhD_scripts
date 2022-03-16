
library(GenomicRanges)
library(rtracklayer)

load("Archive_Year1/BG3_novel_proximal.Rda")
  BG3_proximal_SE <- BG3_novel_proximal[width(BG3_novel_proximal) >= 1000]
load("Archive_Year1/BG3_novel_distal.Rda")
  BG3_distal_SE <- BG3_novel_distal[width(BG3_novel_distal) >= 1000]
load("Archive_Year1/S2_novel_proximal.Rda")
  S2_proximal_SE <- S2_novel_proximal[width(S2_novel_proximal) >= 1000]
load("Archive_Year1/S2_novel_distal.Rda")
  S2_distal_SE <- S2_novel_distal[width(S2_novel_distal) >= 1000]


pdf(file = "/home/jw18713/Project1/Paper_Plots/FigureS2/FigureS2C&D.pdf",
width = 18, height = 8, pointsize = 14)
par(mar=c(5,4.5,5,4), mfrow=c(1,2), cex = 1.2)
par(xpd=NA)

      # BG3 proximal
      hist(log2(width(BG3_proximal_SE)), axes=F, main="BG3 Novel Enhancer Widths",
        xlab="Width (bp)", col="#D55E0080", cex.main=2, cex.lab=1.5, xlim = c(log2(1000),log2(10000)), ylim = c(0,200),
        breaks = 15, ylab=NA)
          axis(1, at=c(log2(1000), log2(2000), log2(4000), log2(6000), log2(8000), log2(10000)),
            labels=c(1000, 2000, 4000, 6000, 8000, 10000))
          axis(2, at=seq(0,200, by = 50), labels = seq(0,200, by = 50))

      #BG3 distal only
      hist(log2(width(BG3_distal_SE)), col="#0072B280", add= T, breaks = 15)

      legend(13.25,150,c("Distal\nOnly", "Proximal"), fill=c("#0072B280", "#D55E0080"), bty="n")

      #S2 proximal
      hist(log2(width(S2_proximal_SE)), axes=F, main="S2 Novel Enhancer Widths",
        xlab="Width (bp)", col="#D55E0080", cex.main=2, cex.lab=1.5, xlim = c(log2(1000),log2(10000)), ylim = c(0,200),
        breaks = 15, ylab=NA)
        axis(1, at=c(log2(1000), log2(2000), log2(4000), log2(6000), log2(8000), log2(10000)),
          labels=c(1000, 2000, 4000, 6000, 8000, 10000))
        axis(2, at=seq(0,200, by = 50), labels = seq(0,200, by = 50))

      #S2 distal only
      hist(log2(width(S2_distal_SE)), col="#0072B280", add = T, breaks = 15)

      legend(13.25,150,c("Distal\nOnly", "Proximal"), fill=c("#0072B280", "#D55E0080"), bty="n")

dev.off()

wilcox.test(width(BG3_proximal_SE), width(BG3_distal_SE))
wilcox.test(width(S2_proximal_SE), width(S2_distal_SE))
