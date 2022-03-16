
library(genomation)
library(GenomicRanges)
library(rtracklayer)

tracks <- read.csv("Review_Paper_Plots/arch_prot_rules/data/STARR_enhancer_widths.csv")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


common <- tracks[tracks$IsCommon == 1, 33:44]
putative <- tracks[tracks$IsCommon == 0, 33:44]

background <- read.csv("Review_Paper_Plots/arch_prot_rules/data/STARR_background.csv")
background <- background[,31:42]


stats_matrix <- matrix(0, 3, ncol(background), dimnames = list(
  c("Putative/Common", "Putative/Background", "Common/Background"),
  c(colnames(background))
))

setwd("~/Review_Paper_Plots/arch_prot_rules/plots/boxplots")
for (i in seq(1, ncol(background))){

  # Boxplot plotting code
  pdf(paste0(colnames(background)[i], "_boxplot_outlilers.pdf"),
    width = 8, height = 6, pointsize = 14)
  par(mar = c(4,5.5,4,8), xpd = T)

  all <- c(putative[,i], common[,i], background[,i])
  limit <- max(abs(range(all)))
  limit <- c(min(all) - limit*0.1, max(all) + limit*0.15)

  boxplot(putative[,i], common[,i], background[,i],
    names = c("Putative", "Common", "Background"),
    #ylim = limit,
    outline = T, col = cbbPalette[2:4], main = colnames(background)[i])
  dev.off()

  # Stats code
  stats_matrix[1,i] <- wilcox.test(putative[,i], common[,i])$p.value
  stats_matrix[2,i] <- wilcox.test(putative[,i], background[,i])$p.value
  stats_matrix[3,i] <- wilcox.test(common[,i], background[,i])$p.value
}

stats_matrix[stats_matrix == 0] <- 2.2e-16

write.csv(stats_matrix, file = "STARR_box_plot_stats_matrix.csv")
