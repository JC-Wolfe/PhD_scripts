
library(GenomicRanges)
library(rtracklayer)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

BG3_contacts <- matrix(0,2,2)
colnames(BG3_contacts) <- c("Distal Only", "Proximal")
rownames(BG3_contacts) <- c("Novel", "Shared")

S2_contacts <- BG3_contacts

BG3_contacts[1,1] <- 9018
BG3_contacts[1,2] <- 8931
BG3_contacts[2,1] <- 1058
BG3_contacts[2,2] <- 1296

S2_contacts[1,1] <- 5508
S2_contacts[1,2] <- 8343
S2_contacts[2,1] <- 799
S2_contacts[2,2] <- 1004

fisher.test(BG3_contacts)
fisher.test(S2_contacts)

BG3_per <- BG3_contacts[,1]
BG3_per[1] <- BG3_contacts[1,1]/(BG3_contacts[1,1]+BG3_contacts[1,2])*100
BG3_per[2] <- BG3_contacts[2,1]/(BG3_contacts[2,1]+BG3_contacts[2,2])*100

S2_per <- S2_contacts[,1]
S2_per[1] <- S2_contacts[1,1]/(S2_contacts[1,1]+S2_contacts[1,2])*100
S2_per[2] <- S2_contacts[2,1]/(S2_contacts[2,1]+S2_contacts[2,2])*100

pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure3/Contact_Percentages.pdf",
width = 14, height = 5, pointsize = 14)
par(mar=c(4,5.5,4,4), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)

  barplot(BG3_per, beside=T, horiz = T, col = cbbPalette[7:6],
      xlim = c(0,100), xlab = "Distal Only Enhancers", xaxt = "n",
      main = "BG3 Contact Distributions")
    axis(1, seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"), las=2)
    lines(c(max(BG3_per) + 2.5, max(BG3_per) + 5), c(1.9,1.9))
    lines(c(max(BG3_per) + 5, max(BG3_per) + 5), c(0.7,1.9))
    lines(c(max(BG3_per) + 2.5, max(BG3_per) + 5), c(0.7,0.7))
    text(x = max(BG3_per) + 7.5, y = 0.7 + (1.9-0.7)/2, "***", srt = 90)

  barplot(S2_per, beside=T, horiz = T, col = cbbPalette[7:6],
      xlim = c(0,100), xlab = "Distal Only Enhancers", xaxt = "n",
      main = "S2 Contact Distributions")
    axis(1, seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"), las=2)
    lines(c(max(S2_per) + 2.5, max(S2_per) + 5), c(1.9,1.9))
    lines(c(max(S2_per) + 5, max(S2_per) + 5), c(0.7,1.9))
    lines(c(max(S2_per) + 2.5, max(S2_per) + 5), c(0.7,0.7))
    text(x = max(S2_per) + 7.5, y = 0.7 + (1.9-0.7)/2, "***", srt = 90)

dev.off()
