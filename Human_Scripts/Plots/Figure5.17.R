library(GenomicRanges)
library(rtracklayer)

setwd("~/Pre_viva/Hi-C/grange_lists")

load("H9_putative_distal_EA.Rda")
load("H9_putative_proximal_EA.Rda")
load("H9_background_distal_EA.Rda")
load("H9_background_proximal_EA.Rda")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7")

setwd("~/Pre_viva/Hi-C/Plots")

hist_density <- function(x){
  h <- hist(x, plot = F, breaks = seq(0, 17, by = 1))
  h$density <- h$counts/sum(h$counts)*100
  return(h)
}


# Ok, so this is my histogram then
pdf(file = "Figure3D.pdf",
width = 14, height = 7, pointsize = 14)
par(mar=c(4.9,3.5,5.5,6), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)

plot(hist_density(log2(H9_putative_distal_EA$max_Expression+1)), freq = F,
  main="H9 Putative Enhancers", xlab="maximum expression of target genes (FPKM)",
  col="#0072B280", cex.main=2, cex.lab=1, axes=F,
  ylim = c(0,100), ylab = NA)
axis(1, at=seq(0,16,by=4), labels=c(0,2^seq(4,16, by = 4)))
axis(2, at=seq(0,100, by = 25), labels = paste0(seq(0, 100, 25), "%"), las=2)

plot(hist_density(log2(H9_putative_proximal_EA$max_Expression+1)), freq = F,
col="#D55E0080", add=T)

#--------------------------Background Expression-------------------------------#
#------------------------------------------------------------------------------#

plot(hist_density(log2(H9_background_distal_EA$max_Expression+1)), freq = F,
  main="H9 Background", xlab="maximum expression of target genes (FPKM)",
  col="#0072B280", cex.main=2, cex.lab=1, axes=F,
  ylim = c(0,100), ylab = NA)
axis(1, at=seq(0,16,by=4), labels=c(0,2^seq(4,16, by = 4)))
axis(2, at=seq(0,100, by = 25), labels = paste0(seq(0, 100, 25), "%"), las=2)

plot(hist_density(log2(H9_background_proximal_EA$max_Expression+1)), freq = F,
col="#D55E0080", add=T)

legend(17,70,c("Distal\nOnly", "Proximal"), fill=c("#0072B280", "#D55E0080"), bty="n")

dev.off()

statsmat <- matrix(0, 1, 2, dimnames = list(c("Proximal/Distal_Only"),
  c("Putative", "Background")))

statsmat[1,1] <- wilcox.test(log2(H9_putative_distal_EA$max_Expression+1),
  log2(H9_putative_proximal_EA$max_Expression+1))$p.value
statsmat[1,2] <- wilcox.test(log2(H9_background_distal_EA$max_Expression+1),
  log2(H9_background_proximal_EA$max_Expression+1))$p.value

write.csv(statsmat, file = "Stats/Figure3D_stats.csv")
