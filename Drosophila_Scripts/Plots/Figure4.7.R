#	Heat map for marks BG3
#	Heat map for marks S2


library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(gplots)
library(lattice)

BG3_results <- get(load("/home/jw18713/Project1/data/model_predictions/BG3_FuzzyGR.Rda")) #Imported as gro
BG3_enhancers <- get(load("~/dm6_promoter_correction/predicted_enhancers/grown_predicted_enhancers_BG3.Rda"))
BG3_starr <- get(load("~/dm6_promoter_correction/starr_seq/starr.Rda"))

S2_results <- get(load("/home/jw18713/Project1/data/model_predictions/S2_FuzzyGR.Rda")) # S2_XAI
S2_enhancers <- get(load("/home/jw18713/dm6_promoter_correction/predicted_enhancers/grown_predicted_enhancers_S2.Rda"))
S2_starr <- get(load("/home/jw18713/dm6_promoter_correction/starr_seq/starr_S2.Rda"))

marks <- mcols(BG3_results)[1:26]

BG3_starr_and_predicted <- reduce(c(BG3_enhancers[BG3_enhancers %over% BG3_starr],
                                BG3_starr[BG3_starr %over% BG3_enhancers]))
BG3_fuzzy <- BG3_enhancers[!BG3_enhancers %over% BG3_starr]
BG3_starr_specific <- BG3_starr[!BG3_starr %over% BG3_enhancers]



BG3_both <- BG3_results[BG3_results %over% BG3_starr_and_predicted]
BG3_fuzzy_only <- BG3_results[BG3_results %over% BG3_fuzzy]
BG3_STARR_only <- BG3_results[BG3_results %over% BG3_starr_specific]
BG3_neither <- BG3_results[!BG3_results %over% BG3_both & !BG3_results %over% BG3_fuzzy_only & !BG3_results %over% BG3_STARR_only]


BG3_map_matrix <- matrix(0, nrow = 4, ncol=length(marks))
colnames(BG3_map_matrix) <- names(marks)
rownames(BG3_map_matrix) <- c("Common", "Putative", "STARR-seq Only", "Neither")

BG3_map_matrix[1,] <- apply(mcols(BG3_both)[1:26], 2, mean)
BG3_map_matrix[2,] <- apply(mcols(BG3_fuzzy_only)[1:26], 2, mean)
BG3_map_matrix[3,] <- apply(mcols(BG3_STARR_only)[1:26], 2, mean)
BG3_map_matrix[4,] <- apply(mcols(BG3_neither)[1:26], 2, mean)

BG3_expected <- apply(mcols(BG3_results)[1:26], 2, mean)

BG3_map_matrix[1,] <- log2(BG3_map_matrix[1,]/BG3_expected)
BG3_map_matrix[2,] <- log2(BG3_map_matrix[2,]/BG3_expected)
BG3_map_matrix[3,] <- log2(BG3_map_matrix[3,]/BG3_expected)
BG3_map_matrix[4,] <- log2(BG3_map_matrix[4,]/BG3_expected)

# Everything below here is S2 stuff

S2_starr_and_predicted <- reduce(c(S2_enhancers[S2_enhancers %over% S2_starr],
                                S2_starr[S2_starr %over% S2_enhancers]))
S2_fuzzy <- S2_enhancers[!S2_enhancers %over% S2_starr]
S2_starr_specific <- S2_starr[!S2_starr %over% S2_enhancers]



S2_both <- S2_results[S2_results %over% S2_starr_and_predicted]
S2_fuzzy_only <- S2_results[S2_results %over% S2_fuzzy]
S2_STARR_only <- S2_results[S2_results %over% S2_starr_specific]
S2_neither <- S2_results[!S2_results %over% S2_both & !S2_results %over% S2_fuzzy_only & !S2_results %over% S2_STARR_only]


S2_map_matrix <- matrix(0, nrow = 4, ncol=length(marks))
colnames(S2_map_matrix) <- names(marks)
rownames(S2_map_matrix) <- c("Common", "Putative", "STARR-seq Only", "Neither")

S2_map_matrix[1,] <- apply(mcols(S2_both)[1:26], 2, mean)
S2_map_matrix[2,] <- apply(mcols(S2_fuzzy_only)[1:26], 2, mean)
S2_map_matrix[3,] <- apply(mcols(S2_STARR_only)[1:26], 2, mean)
S2_map_matrix[4,] <- apply(mcols(S2_neither)[1:26], 2, mean)

S2_expected <- apply(mcols(S2_results)[1:26], 2, mean)

S2_map_matrix[1,] <- log2(S2_map_matrix[1,]/S2_expected)
S2_map_matrix[2,] <- log2(S2_map_matrix[2,]/S2_expected)
S2_map_matrix[3,] <- log2(S2_map_matrix[3,]/S2_expected)
S2_map_matrix[4,] <- log2(S2_map_matrix[4,]/S2_expected)


custom_at <- seq(-0.7, 0.7, by = (1.4/30))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols_contrast <- rev(c(colorRampPalette(c(cbbPalette[7], "white", cbbPalette[6]))(30)))

pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure2/Figure2A_Oct.pdf")
levelplot(t(BG3_map_matrix),
  at = custom_at,
  main = "BG3 Histone Modification Enrichment Change",
  xlab = "log2 Enrichment Difference",
  ylab = "Predicted By",
  scales=list(x=list(rot=90)),
  col.regions = cols_contrast)
dev.off()


pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure2/Figure2B_Oct.pdf")
  levelplot(t(S2_map_matrix),
    at = custom_at,
    main = "S2 Histone Modification Enrichment Change",
    xlab = "log2 Enrichment Difference",
    ylab = "Predicted By",
    scales=list(x=list(rot=90)),
    col.regions = cols_contrast)
dev.off()
