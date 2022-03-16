
library(GenomicRanges)
library(DMRcaller)

load("~/Human_Datasets/H9/Processed_datasets/H9_full.Rda")

H9_putative <- reduce(H9_full[H9_full$classification == "Putative"])
STARR <- reduce(H9_full[H9_full$classification == "STARR"])

d <- distanceToNearest(STARR, H9_putative)
distance <- mcols(d)$distance

pdf(file = "~/Human_Datasets/Plots/Supp/STARR_only_distance_to_nearest.pdf")
hist(log10(distance + 1), main="Distance to Nearest", xlab="Distance",
  col="lightblue", cex.main=2, cex.lab=1.5, axes=F, ylim = c(0, 4000))
axis(1, at=seq(0,6,by=2), labels=c(0,10^seq(2,6, by = 2)))
axis(2, at=seq(0,4000, by = 1000), labels = seq(0,4000, by = 1000))
dev.off()

load("~/Human_Datasets/H9/Processed_datasets/H9_XAI.Rda")
rdat <- GRanges(seqnames = (seqnames(H9_XAI)), ranges = IRanges(start = start(H9_XAI),
end = end(H9_XAI)))
rdat$readsM <- H9_XAI$Conf_Perc_1
rdat$readsN <- 1
rdat$context <- "CG"
rdat$trinucleotide_context <- "CGG"

STARR_avs <- analyseReadsInsideRegionsForCondition(STARR, rdat,
context = "CG", cores = 30)
STARR_conf <- STARR_avs$proportionCG

low_conf <- STARR_avs[STARR_avs$proportionCG <= 0.05]
low_conf <- H9_full[H9_full %over% low_conf]
low_conf$readsM <- H9_XAI$ATA

pdf(file = "~/Human_Datasets/Plots/Supp/STARR_only_avg_confidence.pdf")
hist(STARR_conf, main="STARR-seq Only\nAverage Confidence", xlab="Average Confidence", xlim = c(0, 1),
  col="lightblue", cex.main=2, cex.lab=1.5, axes=F, ylim = c(0, 2000))
axis(1, at=seq(0, 1, by=0.2), labels=seq(0, 1, by = 0.2))
axis(2, at=seq(0,2000, by = 500), labels = seq(0,2000, by = 500))
dev.off()

pdf(file = "~/Human_Datasets/Plots/Supp/STARR_dist_avg_scatter.pdf")
par(mar = c(4, 5.5, 5, 4))
plot(log10(distance + 1), STARR_conf, main = "Average Confidence & Distance", xlab = "log10 Distance", ylab = "Average
  Confidence", axes = F, ylim = c(0,1), xlim = c(0, 6))
axis(1, at=seq(0,6,by=2), labels=c(0,10^seq(2,6, by = 2)))
axis(2, at=seq(0, 1, by=0.2), labels=seq(0, 1, by = 0.2))
dev.off()
