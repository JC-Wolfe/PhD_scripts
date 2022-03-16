library(GenomicRanges)
library(gplots)
library(lattice)

load("~/Human_Datasets/H1/Processed_datasets/H1_full(widths).Rda")
load("~/Human_Datasets/H9/Processed_datasets/H9_full(widths).Rda")

H9_results_matrix <- matrix(0, 5, 6)
colnames(H9_results_matrix) <- c("Intergenic", "Promoter", "Intron", "Exon", "5' utr", "3' utr")
rownames(H9_results_matrix) <- c("Common", "Putative Only", "STARR-seq Only", "Neither", "Whole Genome")

H9_results_matrix[1,1] <- length(which(H9_full$classification == "Common" & H9_full$annotations == "intergenic"))
H9_results_matrix[1,2] <- length(which(H9_full$classification == "Common" & H9_full$annotations == "promoter"))
H9_results_matrix[1,3] <- length(which(H9_full$classification == "Common" & H9_full$annotations == "intron"))
H9_results_matrix[1,4] <- length(which(H9_full$classification == "Common" & H9_full$annotations == "exon"))
H9_results_matrix[1,5] <- length(which(H9_full$classification == "Common" & H9_full$annotations == "five_prime_utr"))
H9_results_matrix[1,6] <- length(which(H9_full$classification == "Common" & H9_full$annotations == "three_prime_utr"))

H9_results_matrix[2,1] <- length(which(H9_full$classification == "Putative" & H9_full$annotations == "intergenic"))
H9_results_matrix[2,2] <- length(which(H9_full$classification == "Putative" & H9_full$annotations == "promoter"))
H9_results_matrix[2,3] <- length(which(H9_full$classification == "Putative" & H9_full$annotations == "intron"))
H9_results_matrix[2,4] <- length(which(H9_full$classification == "Putative" & H9_full$annotations == "exon"))
H9_results_matrix[2,5] <- length(which(H9_full$classification == "Putative" & H9_full$annotations == "five_prime_utr"))
H9_results_matrix[2,6] <- length(which(H9_full$classification == "Putative" & H9_full$annotations == "three_prime_utr"))

H9_results_matrix[3,1] <- length(which(H9_full$classification == "STARR" & H9_full$annotations == "intergenic"))
H9_results_matrix[3,2] <- length(which(H9_full$classification == "STARR" & H9_full$annotations == "promoter"))
H9_results_matrix[3,3] <- length(which(H9_full$classification == "STARR" & H9_full$annotations == "intron"))
H9_results_matrix[3,4] <- length(which(H9_full$classification == "STARR" & H9_full$annotations == "exon"))
H9_results_matrix[3,5] <- length(which(H9_full$classification == "STARR" & H9_full$annotations == "five_prime_utr"))
H9_results_matrix[3,6] <- length(which(H9_full$classification == "STARR" & H9_full$annotations == "three_prime_utr"))

H9_results_matrix[4,1] <- length(which(H9_full$classification == "Neither" & H9_full$annotations == "intergenic"))
H9_results_matrix[4,2] <- length(which(H9_full$classification == "Neither" & H9_full$annotations == "promoter"))
H9_results_matrix[4,3] <- length(which(H9_full$classification == "Neither" & H9_full$annotations == "intron"))
H9_results_matrix[4,4] <- length(which(H9_full$classification == "Neither" & H9_full$annotations == "exon"))
H9_results_matrix[4,5] <- length(which(H9_full$classification == "Neither" & H9_full$annotations == "five_prime_utr"))
H9_results_matrix[4,6] <- length(which(H9_full$classification == "Neither" & H9_full$annotations == "three_prime_utr"))

H9_results_matrix[5,1] <- length(which(H9_full$annotations == "intergenic"))
H9_results_matrix[5,2] <- length(which(H9_full$annotations == "promoter"))
H9_results_matrix[5,3] <- length(which(H9_full$annotations == "intron"))
H9_results_matrix[5,4] <- length(which(H9_full$annotations == "exon"))
H9_results_matrix[5,5] <- length(which(H9_full$annotations == "five_prime_utr"))
H9_results_matrix[5,6] <- length(which(H9_full$annotations == "three_prime_utr"))

H9_hmap <- H9_results_matrix
H9_results_matrix <- t(H9_results_matrix)
H9_results_matrix_aves <- H9_results_matrix

H9_results_matrix_aves[,1] <- H9_results_matrix[,1]/sum(H9_results_matrix[,1])
H9_results_matrix_aves[,2] <- H9_results_matrix[,2]/sum(H9_results_matrix[,2])
H9_results_matrix_aves[,3] <- H9_results_matrix[,3]/sum(H9_results_matrix[,3])
H9_results_matrix_aves[,4] <- H9_results_matrix[,4]/sum(H9_results_matrix[,4])
H9_results_matrix_aves[,5] <- H9_results_matrix[,5]/sum(H9_results_matrix[,5])

H9_hmap_aves <- H9_hmap

H9_hmap_aves[1,] <- H9_hmap[1,]/sum(H9_hmap[1,])
H9_hmap_aves[2,] <- H9_hmap[2,]/sum(H9_hmap[2,])
H9_hmap_aves[3,] <- H9_hmap[3,]/sum(H9_hmap[3,])
H9_hmap_aves[4,] <- H9_hmap[4,]/sum(H9_hmap[4,])
H9_hmap_aves[5,] <- H9_hmap[5,]/sum(H9_hmap[5,])

H9_hmap_comp <- H9_hmap_aves[1:4,]
H9_hmap_comp[1,] <- log2(H9_hmap_comp[1,]/H9_hmap_aves[5,])
H9_hmap_comp[2,] <- log2(H9_hmap_comp[2,]/H9_hmap_aves[5,])
H9_hmap_comp[3,] <- log2(H9_hmap_comp[3,]/H9_hmap_aves[5,])
H9_hmap_comp[4,] <- log2(H9_hmap_comp[4,]/H9_hmap_aves[5,])

H1_results_matrix <- matrix(0, 5, 6)
colnames(H1_results_matrix) <- c("Intergenic", "Promoter", "Intron", "Exon", "5' utr", "3' utr")
rownames(H1_results_matrix) <- c("Common", "H1 Specific", "H9 Specific", "Neither", "Whole Genome")

H1_results_matrix[1,1] <- length(which(H1_full$classification == "Common" & H1_full$annotations == "intergenic"))
H1_results_matrix[1,2] <- length(which(H1_full$classification == "Common" & H1_full$annotations == "promoter"))
H1_results_matrix[1,3] <- length(which(H1_full$classification == "Common" & H1_full$annotations == "intron"))
H1_results_matrix[1,4] <- length(which(H1_full$classification == "Common" & H1_full$annotations == "exon"))
H1_results_matrix[1,5] <- length(which(H1_full$classification == "Common" & H1_full$annotations == "five_prime_utr"))
H1_results_matrix[1,6] <- length(which(H1_full$classification == "Common" & H1_full$annotations == "three_prime_utr"))

H1_results_matrix[2,1] <- length(which(H1_full$classification == "H1 Specific" & H1_full$annotations == "intergenic"))
H1_results_matrix[2,2] <- length(which(H1_full$classification == "H1 Specific" & H1_full$annotations == "promoter"))
H1_results_matrix[2,3] <- length(which(H1_full$classification == "H1 Specific" & H1_full$annotations == "intron"))
H1_results_matrix[2,4] <- length(which(H1_full$classification == "H1 Specific" & H1_full$annotations == "exon"))
H1_results_matrix[2,5] <- length(which(H1_full$classification == "H1 Specific" & H1_full$annotations == "five_prime_utr"))
H1_results_matrix[2,6] <- length(which(H1_full$classification == "H1 Specific" & H1_full$annotations == "three_prime_utr"))

H1_results_matrix[3,1] <- length(which(H1_full$classification == "H9 Specific" & H1_full$annotations == "intergenic"))
H1_results_matrix[3,2] <- length(which(H1_full$classification == "H9 Specific" & H1_full$annotations == "promoter"))
H1_results_matrix[3,3] <- length(which(H1_full$classification == "H9 Specific" & H1_full$annotations == "intron"))
H1_results_matrix[3,4] <- length(which(H1_full$classification == "H9 Specific" & H1_full$annotations == "exon"))
H1_results_matrix[3,5] <- length(which(H1_full$classification == "H9 Specific" & H1_full$annotations == "five_prime_utr"))
H1_results_matrix[3,6] <- length(which(H1_full$classification == "H9 Specific" & H1_full$annotations == "three_prime_utr"))

H1_results_matrix[4,1] <- length(which(H1_full$classification == "Neither" & H1_full$annotations == "intergenic"))
H1_results_matrix[4,2] <- length(which(H1_full$classification == "Neither" & H1_full$annotations == "promoter"))
H1_results_matrix[4,3] <- length(which(H1_full$classification == "Neither" & H1_full$annotations == "intron"))
H1_results_matrix[4,4] <- length(which(H1_full$classification == "Neither" & H1_full$annotations == "exon"))
H1_results_matrix[4,5] <- length(which(H1_full$classification == "Neither" & H1_full$annotations == "five_prime_utr"))
H1_results_matrix[4,6] <- length(which(H1_full$classification == "Neither" & H1_full$annotations == "three_prime_utr"))

H1_results_matrix[5,1] <- length(which(H1_full$annotations == "intergenic"))
H1_results_matrix[5,2] <- length(which(H1_full$annotations == "promoter"))
H1_results_matrix[5,3] <- length(which(H1_full$annotations == "intron"))
H1_results_matrix[5,4] <- length(which(H1_full$annotations == "exon"))
H1_results_matrix[5,5] <- length(which(H1_full$annotations == "five_prime_utr"))
H1_results_matrix[5,6] <- length(which(H1_full$annotations == "three_prime_utr"))

H1_hmap <- H1_results_matrix
H1_results_matrix <- t(H1_results_matrix)
H1_results_matrix_aves <- H1_results_matrix

H1_results_matrix_aves[,1] <- H1_results_matrix[,1]/sum(H1_results_matrix[,1])
H1_results_matrix_aves[,2] <- H1_results_matrix[,2]/sum(H1_results_matrix[,2])
H1_results_matrix_aves[,3] <- H1_results_matrix[,3]/sum(H1_results_matrix[,3])
H1_results_matrix_aves[,4] <- H1_results_matrix[,4]/sum(H1_results_matrix[,4])
H1_results_matrix_aves[,5] <- H1_results_matrix[,5]/sum(H1_results_matrix[,5])

H1_hmap_aves <- H1_hmap

H1_hmap_aves[1,] <- H1_hmap[1,]/sum(H1_hmap[1,])
H1_hmap_aves[2,] <- H1_hmap[2,]/sum(H1_hmap[2,])
H1_hmap_aves[3,] <- H1_hmap[3,]/sum(H1_hmap[3,])
H1_hmap_aves[4,] <- H1_hmap[4,]/sum(H1_hmap[4,])
H1_hmap_aves[5,] <- H1_hmap[5,]/sum(H1_hmap[5,])

H1_hmap_comp <- H1_hmap_aves[1:4,]
H1_hmap_comp[1,] <- log2(H1_hmap_comp[1,]/H1_hmap_aves[5,])
H1_hmap_comp[2,] <- log2(H1_hmap_comp[2,]/H1_hmap_aves[5,])
H1_hmap_comp[3,] <- log2(H1_hmap_comp[3,]/H1_hmap_aves[5,])
H1_hmap_comp[4,] <- log2(H1_hmap_comp[4,]/H1_hmap_aves[5,])


cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
custom_at <- seq(-3, 3, by = (6/60))
cols_contrast <- rev(c(colorRampPalette(c(cbbPalette[6], "white", cbbPalette[5]))(60)))

H9_hmap_comp[H9_hmap_comp < -3] <- -3
H9_hmap_comp[H9_hmap_comp > 3] <- 3

H1_hmap_comp[H1_hmap_comp < -3] <- -3
H1_hmap_comp[H1_hmap_comp > 3] <- 3

pdf("/home/jw18713/Human_Datasets/Plots/Figure2/H1&H9_Annotations(width).pdf",
width = 14, height = 8, pointsize = 14)
par(mar=c(6,4,5,7), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)
x <- barplot(H9_results_matrix_aves, col=cbbPalette,
  ylab="Percentage of Results in Region", main="H9 & STARR-seq\nAnnotation Comparisons",
  yaxt = "none", xaxt = "none")
labs <- colnames(H9_results_matrix_aves)
text(cex=1, x=x, y=-0.04, labs, adj = 1, srt=45)
axis(2, seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10), "%"), las=2)
legend(6.25,0.75,rownames(H9_results_matrix_aves), fill=cbbPalette, bty = "n")

x <- barplot(H1_results_matrix_aves, col=cbbPalette,
  ylab="Percentage of Results in Region", main="H1 & H9\nAnnotation Comparisons",
  yaxt = "none", xaxt = "none")
labs <- colnames(H1_results_matrix_aves)
text(cex=1, x=x, y=-0.04, labs, adj = 1, srt=45)
axis(2, seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10), "%"), las=2)
axis(2, seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10), "%"), las=2)
legend(6.25,0.75,rownames(H1_results_matrix_aves), fill=cbbPalette, bty = "n")
dev.off()


pdf(file = "/home/jw18713/Human_Datasets/Plots/Figure2/H9_Annotation_Heatmap(width).pdf")
  levelplot(t(H9_hmap_comp),
    at = custom_at,
    main = "Predicted Region\nAnnotation Comparisons H9",
    xlab = "log2 Enrichment Difference",
    ylab = "Classification",
    scales=list(x=list(rot=90)),
    col.regions = cols_contrast)
dev.off()

pdf(file = "/home/jw18713/Human_Datasets/Plots/Figure2/H1_Annotation_Heatmap(width).pdf")
  levelplot(t(H1_hmap_comp),
    at = custom_at,
    main = "H1 & H9 Predicted Region\nAnnotation Comparisons",
    xlab = "log2 Enrichment Difference",
    ylab = "Classification",
    scales=list(x=list(rot=90)),
    col.regions = cols_contrast)
dev.off()
