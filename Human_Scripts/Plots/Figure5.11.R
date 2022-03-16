library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(gplots)
library(lattice)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
load("~/Human_Datasets/H1/Processed_datasets/H1_full(widths).Rda")
load("~/Human_Datasets/H9/Processed_datasets/H9_full(widths).Rda")
H1_full$ATAC_signalValues <- NULL

zero_to_NA <- function(x){
  x[x == 0] <- NA
  return(x)
}

NA_cols <- apply(mcols(H9_full)[1:27], 2, zero_to_NA)
mcols(H9_full)[1:27] <- NA_cols

NA_cols <- apply(mcols(H1_full)[1:27], 2, zero_to_NA)
mcols(H1_full)[1:27] <- NA_cols

H9_common <- H9_full[H9_full$classification == "Common"]
H9_putative <- H9_full[H9_full$classification == "Putative"]
H9_neither <- H9_full[H9_full$classification == "Neither"]
H9_starr <- H9_full[H9_full$classification == "STARR"]

H1_common <- H1_full[H1_full$classification == "Common"]
H1_H1_specific <- H1_full[H1_full$classification == "H1 Specific"]
H1_H9_specific <- H9_full[H1_full$classification == "H9 Specific"]
H1_neither <- H1_full[H9_full$classification == "Neither"]

H9_map_matrix <- matrix(0, nrow = 4, ncol=31)
colnames(H9_map_matrix) <- names(mcols(H9_full)[1:31])
rownames(H9_map_matrix) <- c("Common", "Putative", "STARR-seq", "Neither")

H9_map_matrix[1,] <- apply(mcols(H9_common)[1:31], 2, mean, na.rm = T)
H9_map_matrix[2,] <- apply(mcols(H9_putative)[1:31], 2, mean, na.rm = T)
H9_map_matrix[3,] <- apply(mcols(H9_starr)[1:31], 2, mean, na.rm = T)
H9_map_matrix[4,] <- apply(mcols(H9_neither)[1:31], 2, mean, na.rm = T)

H9_expected <- apply(mcols(H9_full)[1:31], 2, mean, na.rm = T)

H9_map_matrix[1,] <- log2(H9_map_matrix[1,]/H9_expected)
H9_map_matrix[2,] <- log2(H9_map_matrix[2,]/H9_expected)
H9_map_matrix[3,] <- log2(H9_map_matrix[3,]/H9_expected)
H9_map_matrix[4,] <- log2(H9_map_matrix[4,]/H9_expected)

# H1 stuff
H1_map_matrix <- matrix(0, nrow = 4, ncol=31)
colnames(H1_map_matrix) <- names(mcols(H1_full)[1:31])
rownames(H1_map_matrix) <- c("Common", "H1 Specific", "H9 Specific", "Neither")

H1_map_matrix[1,] <- apply(mcols(H1_common)[1:31], 2, mean, na.rm = T)
H1_map_matrix[2,] <- apply(mcols(H1_H1_specific)[1:31], 2, mean, na.rm = T)
H1_map_matrix[3,] <- apply(mcols(H1_H9_specific)[1:31], 2, mean, na.rm = T)
H1_map_matrix[4,] <- apply(mcols(H1_neither)[1:31], 2, mean, na.rm = T)

H1_expected <- apply(mcols(H1_full)[1:31], 2, mean, na.rm = T)

H1_map_matrix[1,] <- log2(H1_map_matrix[1,]/H1_expected)
H1_map_matrix[2,] <- log2(H1_map_matrix[2,]/H1_expected)
H1_map_matrix[3,] <- log2(H1_map_matrix[3,]/H9_expected)
H1_map_matrix[4,] <- log2(H1_map_matrix[4,]/H1_expected)

custom_at <- seq(-1.5, 1.5, by = (3/30))

H9_map_matrix[H9_map_matrix < -1.5] <- -1.5
H9_map_matrix[H9_map_matrix > 1.5] <- 1.5

H1_map_matrix[H1_map_matrix < -1.5] <- -1.5
H1_map_matrix[H1_map_matrix > 1.5] <- 1.5

cols_contrast <- rev(c(colorRampPalette(c(cbbPalette[7], "white", cbbPalette[6]))(60)))

pdf(file = "~/Human_Datasets/Plots/Figure2/H9_Heatmap(widths).pdf")
levelplot(t(H9_map_matrix),
  at = custom_at,
  main = "H9 Histone Modification Enrichment Change",
  xlab = "log2 Enrichment Difference",
  ylab = "Classification",
  scales=list(x=list(rot=90)),
  col.regions = cols_contrast)
dev.off()

pdf(file = "~/Human_Datasets/Plots/Figure2/H1vsH9_Heatmap(widths).pdf")
levelplot(t(H1_map_matrix),
  at = custom_at,
  main = "H1 vs H9 Histone Modification Enrichment Change",
  xlab = "log2 Enrichment Difference",
  ylab = "Classification",
  scales=list(x=list(rot=90)),
  col.regions = cols_contrast)
dev.off()
