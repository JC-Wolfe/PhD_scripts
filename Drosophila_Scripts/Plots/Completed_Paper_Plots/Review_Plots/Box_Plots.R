
library(genomation)
library(GenomicRanges)
library(rtracklayer)

load("Review_Paper_Plots/rda_stuff/widths.Rda")
load("Review_Paper_Plots/rda_stuff/new_order.Rda")
scaled_200 <- read.csv("Review_Paper_Plots/arch_prot_rules/data/ChIP_enrichment_95perc_quartile_table_enhancer_width_window_ranked.csv")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


common <- scaled_200[scaled_200$IsCommon == 1,-c(1,15:22,25)]
putative <- scaled_200[scaled_200$IsCommon == 0, -c(1,15:22,25)]

background <- read.csv("Review_Paper_Plots/arch_prot_rules/data/ChIP_enrichment_95perc_quartile_table_2Kb_genome_ranked.csv")
background <- background[,-c(1:2,16:23)]

boxplot(putative$Nipped.B, common$Nipped.B, background$Nipped.B,
  outline = F, col = cbbPalette[2:4])

stats_matrix <- matrix(0, 3, ncol(background), dimnames = list(
  c("Putative/Common", "Putative/Background", "Common/Background"),
  c(colnames(background))
))

setwd("~/Review_Paper_Plots/arch_prot_rules/plots/boxplots")
for (i in seq(1, ncol(background))){

  # Boxplot plotting code
  pdf(paste0(colnames(background)[i], "_boxplot.pdf"),
    width = 8, height = 6, pointsize = 14)
  par(mar = c(4,5.5,4,8), xpd = T)
  boxplot(putative[,i], common[,i], background[,i],
    outline = F, col = cbbPalette[2:4], main = colnames(background)[i],
    ylim = c(0,1))
  legend(3.75, 0.7, legend = c("Putative", "Common", "Background"),
         fill = cbbPalette[2:4], title = "Classification", bty="n")
  dev.off()

  # Stats code
  stats_matrix[1,i] <- wilcox.test(putative[,i], common[,i])$p.value
  stats_matrix[2,i] <- wilcox.test(putative[,i], background[,i])$p.value
  stats_matrix[3,i] <- wilcox.test(common[,i], background[,i])$p.value
}

stats_matrix[stats_matrix == 0] <- 2.2e-16

write.csv(stats_matrix, file = "box_plot_stats_matrix.csv")
