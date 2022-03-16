library(GenomicRanges)
library(rtracklayer)

setwd("~/Pre_viva/Hi-C")
# Loading in contact lists
load("contact_Rda/H9_putative.Rda")
load("contact_Rda/H9_common.Rda")

# Getting widths of distal only and proximal enhancers
# Proximal is proximal only [[3]] & distal and proximal [[1]]
H9_putative_DO <- width(H9_putative[[2]])
H9_putative_PR <- width(c(H9_putative[[1]],H9_putative[[3]]))

# Getting widths of common enhancers (same combinations and indices)
H9_common_DO <- width(H9_common[[2]])
H9_common_PR <- width(c(H9_common[[1]],H9_common[[3]]))

hist_density <- function(x, n = 20){
  h <- hist(x, plot = F, breaks = seq(2, 5, by = (5-2) / 12))
  h$density <- h$counts/sum(h$counts)*100
  return(h)
}

setwd("~/Pre_viva/Hi-C/Plots")
# Plotting Figure 3B
pdf(file = "Figure3B_1.pdf",
width = 9, height = 8, pointsize = 14)
par(mar=c(5,5.5,5,6), cex = 1.2)
par(xpd = NA)

plot(hist_density(log10(H9_putative_PR)), freq = F, xlim = c(2,5),
  ylim = c(0,25), axes=F, main="H9 Putative Enhancer Widths", xlab="Width (bp)",
  col="#D55E0080", cex.main=2, cex.lab=1.5, ylab="Density")
axis(1, at=c(log10(100), log10(1000), log10(10000), log10(100000)),
  labels=c("100", "1000", "10000", "100000"))
axis(2, seq(0, 25, 5), labels = paste0(seq(0, 25, 5), "%"), las=2)

plot(hist_density(log10(H9_putative_DO)), freq = F, col="#0072B280", add= T)

legend(5,17.5,c("Distal\nOnly", "Proximal"), fill=c("#0072B280", "#D55E0080"), bty="n")

dev.off()


# Plotting Figure 3B
pdf(file = "Figure3B_2.pdf",
width = 9, height = 8, pointsize = 14)
par(mar=c(5,5.5,5,6), cex = 1.2)
par(xpd = NA)

plot(hist_density(log10(H9_common_PR)), freq = F, xlim = c(2,5),
  ylim = c(0,75), axes=F, main="H9 Common Enhancer Widths", xlab="Width (bp)",
  col="#D55E0080", cex.main=2, cex.lab=1.5, ylab="Density")
axis(1, at=c(log10(100), log10(1000), log10(10000), log10(100000)),
  labels=c("100", "1000", "10000", "100000"))
axis(2, seq(0, 75, 15), labels = paste0(seq(0, 75, 15), "%"), las=2)

plot(hist_density(log10(H9_common_DO)), freq = F, col="#0072B280", add= T)

legend(5,40,c("Distal\nOnly", "Proximal"), fill=c("#0072B280", "#D55E0080"), bty="n")

dev.off()

statsmat <- matrix(0, 1, 2, dimnames = list(c("Proximal/Distal_Only"),
  c("Putative", "Common")))

statsmat[1,1] <- wilcox.test(H9_putative_PR, H9_putative_DO)$p.value
statsmat[1,2] <- wilcox.test(H9_common_PR, H9_common_DO)$p.value

write.csv(statsmat, file = "Stats/Figure3B_stats.csv")




write.csv(statsmat, file = "Stats/Figure3B_stats.csv")
