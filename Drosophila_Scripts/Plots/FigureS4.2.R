library(GenomicRanges)
library(rtracklayer)

load("~/Paper_replots/Objects/Expression_lists/S2_starr_proximal.Rda")
load("~/Paper_replots/Objects/Expression_lists/S2_starr_distal_only.Rda")
load("~/Paper_replots/Objects/Expression_lists/BG3_starr_proximal.Rda")
load("~/Paper_replots/Objects/Expression_lists/BG3_starr_distal_only.Rda")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7")

# Ok, so this is my histogram then
pdf(file = "~/Paper_replots/Plots/Supplementary/Figure_S5.pdf",
width = 20, height = 8, pointsize = 14)
par(mar=c(4.9,5.5,5.5,9), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)

hist(log2(BG3_starr_distal_only$max_Expression+1),
  main="BG3 STARR-seq only enhancers", xlab="maximum expression of target genes (FPKM)",
  col="#0072B280", cex.main=2, cex.lab=1, axes=F,
  ylim = c(0,300), breaks = 16, ylab = NA)
axis(1, at=seq(0,16,by=4), labels=c(0,2^seq(4,16, by = 4)))
axis(2, at=seq(0,300, by = 100), labels = seq(0,300, by = 100), las=2)

hist(log2(BG3_starr_proximal$max_Expression+1), col="#D55E0080",
add=T, breaks = 16)

# S2 stuff

hist(log2(S2_starr_distal_only$max_Expression+1),
  main="S2 STARR-seq only enhancers",
  xlab="maximum expression of target genes (FPKM)", col="#0072B280", cex.main=2, cex.lab=1, axes=F,
  ylim = c(0,300), breaks = 16, ylab = NA)
axis(1, at=seq(0,16,by=4), labels=c(0,2^seq(4,16, by = 4)))
axis(2, at=seq(0,300, by = 100), labels = seq(0,300, by = 100), las=2)

hist(log2(S2_starr_proximal$max_Expression+1), col="#D55E0080",
  add=T, breaks = 16)

legend(16,200,c("Distal only", "Proximal"), fill=c("#0072B280", "#D55E0080"),
  bty='n')

dev.off()

statsmatrix_S5 <- matrix(0, 1, 2, dimnames = list(c("STARR-seq Expressions"), c("BG3", "S2")))

# Stats for figure legends
statsmatrix_S5[1,1] <- wilcox.test(log2(BG3_starr_proximal$max_Expression+1),
  log2(BG3_starr_distal_only$max_Expression+1))$p.value

statsmatrix_S5[1,2] <- wilcox.test(log2(S2_starr_proximal$max_Expression+1),
  log2(S2_starr_distal_only$max_Expression+1))$p.value

write.csv(statsmatrix_S5, file = "~/Paper_replots/Stats/S5_stats.csv")
