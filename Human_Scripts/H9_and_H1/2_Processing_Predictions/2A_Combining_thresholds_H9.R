library(GenomicRanges)
library(DMRcaller)

load("~/Human_Datasets/H9/Processed_datasets/H9_XAI.Rda")
over_thresh <- H9_XAI[H9_XAI$Conf_Perc_1 >= 0.8]

rdat <- GRanges(seqnames = (seqnames(H9_XAI)), ranges = IRanges(start = start(H9_XAI),
end = end(H9_XAI)))
rdat$readsM <- H9_XAI$Conf_Perc_1
rdat$readsN <- 1
rdat$context <- "CG"
rdat$trinucleotide_context <- "CGG"

unreduced <- over_thresh
red_500 <- reduce(over_thresh, min.gapwidth = 500)
red_1000 <- reduce(over_thresh, min.gapwidth = 1000)
red_1500 <- reduce(over_thresh, min.gapwidth = 1500)
red_2000 <- reduce(over_thresh, min.gapwidth = 2000)
red_2500 <- reduce(over_thresh, min.gapwidth = 2500)
red_3000 <- reduce(over_thresh, min.gapwidth = 3000)

counts <- c(length(unreduced), length(red_500), length(red_1000),
length(red_1500), length(red_2000), length(red_2500), length(red_3000))

unreduced_avs <- analyseReadsInsideRegionsForCondition(unreduced, rdat,
  context = "CG", cores = 30)
mean_ur <- mean(unreduced_avs$proportionCG)
width_ur <- mean(width(unreduced))

red_500_avs <- analyseReadsInsideRegionsForCondition(red_500, rdat,
  context = "CG", cores = 30)
mean_500 <- mean(red_500_avs$proportionCG)
width_500 <- mean(width(red_500))

red_1000_avs <- analyseReadsInsideRegionsForCondition(red_1000, rdat,
  context = "CG", cores = 30)
mean_1000 <- mean(red_1000_avs$proportionCG)
width_1000 <- mean(width(red_1000))

red_1500_avs <- analyseReadsInsideRegionsForCondition(red_1500, rdat,
  context = "CG", cores = 30)
mean_1500 <- mean(red_1500_avs$proportionCG)
width_1500 <- mean(width(red_1500))

red_2000_avs <- analyseReadsInsideRegionsForCondition(red_2000, rdat,
  context = "CG", cores = 30)
mean_2000 <- mean(red_2000_avs$proportionCG)
width_2000 <- mean(width(red_2000))

red_2500_avs <- analyseReadsInsideRegionsForCondition(red_2500, rdat,
  context = "CG", cores = 30)
mean_2500 <- mean(red_2500_avs$proportionCG)
width_2500 <- mean(width(red_2500))

red_3000_avs <- analyseReadsInsideRegionsForCondition(red_3000, rdat,
  context = "CG", cores = 30)
mean_3000 <- mean(red_3000_avs$proportionCG)
width_3000 <- mean(width(red_3000))

averages <- c(mean_ur, mean_500, mean_1000, mean_1500, mean_2000, mean_2500,
  mean_3000)

widths <- c(width_ur, width_500, width_1000, width_1500, width_2000,
  width_2500, width_3000)

H9_reduction_matrix <- matrix(c(counts, averages, widths), 3, 7, byrow = T,
  dimnames = list(c("Count", "Mean_Conf_Perc_1", "Mean Width"),
  c("0", "500", "1000", "1500", "2000", "2500", "3000")))

save(H9_reduction_matrix,
  file = "~/Human_Datasets/H9/Processed_datasets/H9_reduction_matrix.Rda")
