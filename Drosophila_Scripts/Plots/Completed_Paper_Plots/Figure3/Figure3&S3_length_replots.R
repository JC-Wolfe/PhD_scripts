library(GenomicRanges)
library(rtracklayer)

# Loading in contact lists
load("~/dm6_promoter_correction/Class_split_glists/BG3_putative_list.Rda")
load("~/dm6_promoter_correction/Class_split_glists/BG3_shared_list.Rda")
load("~/dm6_promoter_correction/Class_split_glists/S2_putative_list.Rda")
load("~/dm6_promoter_correction/Class_split_glists/S2_shared_list.Rda")

# Getting widths of distal only and proximal enhancers
# Proximal is proximal only [[3]] & distal and proximal [[1]]
BG3_putative_DO <- width(BG3_putative[[2]])
BG3_putative_PR <- width(c(BG3_putative[[1]],BG3_putative[[3]]))
S2_putative_DO <- width(S2_putative[[2]])
S2_putative_PR <- width(c(S2_putative[[1]],S2_putative[[3]]))

# Getting widths of common enhancers (same combinations and indices)
BG3_common_DO <- width(BG3_shared[[2]])
BG3_common_PR <- width(c(BG3_shared[[1]],BG3_shared[[3]]))
S2_common_DO <- width(S2_shared[[2]])
S2_common_PR <- width(c(S2_shared[[1]],S2_shared[[3]]))

# Plotting Figure 3B
pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure3/Figure3B_redone.pdf",
width = 18, height = 8, pointsize = 14)
par(mar=c(5,5.5,5,6), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)

      # BG3 proximal
      hist(log2(BG3_putative_PR), axes=F, main="BG3 Novel Enhancer Widths",
        xlab="Width (bp)", col="#D55E0080", cex.main=2, cex.lab=1.5, xlim = c(3,14), ylim = c(0,1750),
        breaks = 20, ylab=NA)
          axis(1, at=c(log2(10), log2(50), log2(150), log2(400), log2(1000), log2(4000), log2(16000)),
              labels=c(10, 50, 150, 400, 1000, 4000, 16000))
          axis(2, at=seq(0,1750, by = 250), labels = seq(0,1750, by = 250), las=2)
          lines(x=c(log2(50),log2(50)), y=c(0,1750), col="#D55E00", lty=2, lwd=2)
          lines(x=c(log2(1000),log2(1000)), y=c(0,1750), col="#D55E00", lty=2, lwd=2)
          text(4.25, 1700, "Fragments")
          text((log2(1000) - log2(50))/2 + log2(50), 1700, "Enhancers")
          text((14 - log2(1000))/2 + log2(1000), 1700, "Super\nEnhancers")

      #BG3 distal only
      hist(log2(BG3_putative_DO), col="#0072B280", add= T, breaks = 20)

      legend(14.5,1200,c("Distal\nOnly", "Proximal"), fill=c("#0072B280", "#D55E0080"), bty="n")

      #S2 proximal
      hist(log2(S2_putative_PR), axes=F, main="S2 Novel Enhancer Widths",
        xlab="Width (bp)", col="#D55E0080", cex.main=2, cex.lab=1.5, xlim = c(3,14), ylim = c(0,1750),
        breaks = 20, ylab=NA)
          axis(1, at=c(log2(10), log2(50), log2(150), log2(400), log2(1000), log2(4000), log2(16000)),
              labels=c(10, 50, 150, 400, 1000, 4000, 16000))
          axis(2, at=seq(0,1750, by = 250), labels = seq(0,1750, by = 250), las=2)
          lines(x=c(log2(50),log2(50)), y=c(0,1750), col="#D55E00", lty=2, lwd=2)
          lines(x=c(log2(1000),log2(1000)), y=c(0,1750), col="#D55E00", lty=2, lwd=2)
          text(4.25, 1700, "Fragments")
          text((log2(1000) - log2(50))/2 + log2(50), 1700, "Enhancers")
          text((14 - log2(1000))/2 + log2(1000), 1700, "Super\nEnhancers")

      #S2 distal only
      hist(log2(S2_putative_DO), col="#0072B280", add = T, breaks = 20)
      legend(14.5,1200,c("Distal\nOnly", "Proximal"), fill=c("#0072B280", "#D55E0080"), bty="n")

dev.off()

# Stats for figure legends
wilcox.test(BG3_putative_PR, BG3_putative_DO)
wilcox.test(S2_putative_PR, S2_putative_DO)

# Plotting Figure S3B
pdf(file = "/home/jw18713/Project1/Paper_Plots/FigureS3/FigureS3B_redone.pdf",
width = 18, height = 8, pointsize = 14)
par(mar=c(5,5.5,5,6), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)

      # BG3 proximal
      hist(log2(BG3_common_PR), axes=F, main="BG3 Novel Enhancer Widths",
        xlab="Width (bp)", col="#D55E0080", cex.main=2, cex.lab=1.5, xlim = c(3,14), ylim = c(0,300),
        breaks = 20, ylab=NA)
          axis(1, at=c(log2(10), log2(50), log2(150), log2(400), log2(1000), log2(4000), log2(16000)),
              labels=c(10, 50, 150, 400, 1000, 4000, 16000))
          axis(2, at=seq(0,300, by = 100), labels = seq(0,300, by = 100), las=2)
          lines(x=c(log2(50),log2(50)), y=c(0,300), col="#D55E00", lty=2, lwd=2)
          lines(x=c(log2(1000),log2(1000)), y=c(0,300), col="#D55E00", lty=2, lwd=2)
          text(4.25, 290, "Fragments")
          text((log2(1000) - log2(50))/2 + log2(50), 290, "Enhancers")
          text((14 - log2(1000))/2 + log2(1000), 290, "Super\nEnhancers")

      #BG3 distal only
      hist(log2(BG3_common_DO), col="#0072B280", add= T, breaks = 20)

      legend(14.5,200,c("Distal\nOnly", "Proximal"), fill=c("#0072B280", "#D55E0080"), bty="n")

      #S2 proximal
      hist(log2(S2_common_PR), axes=F, main="S2 Novel Enhancer Widths",
        xlab="Width (bp)", col="#D55E0080", cex.main=2, cex.lab=1.5, xlim = c(3,14), ylim = c(0,300),
        breaks = 20, ylab=NA)
          axis(1, at=c(log2(10), log2(50), log2(150), log2(400), log2(1000), log2(4000), log2(16000)),
              labels=c(10, 50, 150, 400, 1000, 4000, 16000))
          axis(2, at=seq(0,300, by = 100), labels = seq(0,300, by = 100), las=2)
          lines(x=c(log2(50),log2(50)), y=c(0,300), col="#D55E00", lty=2, lwd=2)
          lines(x=c(log2(1000),log2(1000)), y=c(0,300), col="#D55E00", lty=2, lwd=2)
          text(4.25, 290, "Fragments")
          text((log2(1000) - log2(50))/2 + log2(50), 290, "Enhancers")
          text((14 - log2(1000))/2 + log2(1000), 290, "Super\nEnhancers")

      #S2 distal only
      hist(log2(S2_common_DO), col="#0072B280", add = T, breaks = 20)
      legend(14.5,200,c("Distal\nOnly", "Proximal"), fill=c("#0072B280", "#D55E0080"), bty="n")

dev.off()

# Stats for figure legends
wilcox.test(BG3_common_PR, BG3_common_DO)
wilcox.test(S2_common_PR, S2_common_DO)
