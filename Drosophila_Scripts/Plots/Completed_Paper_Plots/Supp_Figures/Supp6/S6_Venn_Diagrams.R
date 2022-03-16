
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(VennDiagram)
library(seqLogo)
library(GenomicRanges)
library(BSgenome.Dmelanogaster.UCSC.dm6)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
library(PWMEnrich)
library(PWMEnrich.Dmelanogaster.background)
library(MotifDb)
library(gridExtra)
# Load the required packages
require(ggplot2)
require(ggseqlogo)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


BG3_put <- read.table("/home/jw18713/Project1/Paper_Tables/BG3_putative_diff_enriched_TFs.dat", stringsAsFactors = F)[,1]
S2_put <- read.table("/home/jw18713/Project1/Paper_Tables/S2_putative_diff_enriched_TFs.dat", stringsAsFactors = F)[,1]
BG3_common <- read.table("/home/jw18713/Project1/Paper_Tables/BG3_common_diff_enriched_TFs.dat", stringsAsFactors = F)[,1]
S2_common <- read.table("/home/jw18713/Project1/Paper_Tables/S2_common_diff_enriched_TFs.dat", stringsAsFactors = F)[,1]

put_all <- unique(c(BG3_put, S2_put))
common_all <- unique(c(BG3_common, S2_common))
shared <- put_all[put_all%in%common_all]

put_both <- unique(BG3_put[BG3_put%in%S2_put])
common_both <- unique(BG3_common[BG3_common%in%S2_common])
shared2 <- put_both[put_both%in%common_both]

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

generatePairwiseVenn <- function(set1, set2, categories, cols=c("#0072B2", "#D55E00"), cat.pos = c(-20, 20), cat.dist = c(0.07, 0.07), human = TRUE){

  v <- draw.pairwise.venn(area1 = length(set1), area2 = length(set2), cross.area = sum(set1%in%set2),
                          category = categories, col = "transparent", fill = cols,
                          alpha = 0.6, label.col = rep("black", 3), cex = 1.2,
                          cat.col =  cols, cat.cex = 1.4,
                          cat.pos = cat.pos,
                          cat.dist = cat.dist,
                          margin = 0.2,
                          euler.d =TRUE, scaled = T
  )
  if(human){
    for(i in 5:7){
      v[[i]]$label  <- as.vector(sciNotation(as.numeric(v[[i]]$label)))
    }
  }

  return(v)
}



v1 <- generatePairwiseVenn(BG3_put, S2_put, c("BG3", "S2"),
                           cols=c(cbbPalette[6], cbbPalette[7]),
                           cat.pos = c(-30, 150),
                           cat.dist = c(0.075, 0.075), human = F)


v2 <- generatePairwiseVenn(BG3_common, S2_common, c("BG3", "S2"),
                           cols=c(cbbPalette[6], cbbPalette[7]),
                           cat.pos = c(-30, 150),
                           cat.dist = c(0.075, 0.075), human = F)


pdf("/home/jw18713/Project1/Paper_Plots/FigureS6/FigureS6_venn.pdf", width=8, height=5,pointsize = 10);
par(mar=c(0, 0, 0, 0)+0.1)
pushViewport(plotViewport(layout=grid.layout(1, 2),gp=gpar(cex=1.0)))
pushViewport(plotViewport(layout.pos.row=1, layout.pos.col=1))
grid.draw(v1)
grid.text(bquote(bold(""~"Putative")), y=0.93, gp = gpar(cex=1.2))
grid.text(LETTERS[1], y=0.93, x=0.07, gp = gpar(cex=1.6))
popViewport()
pushViewport(plotViewport(layout.pos.row=1, layout.pos.col=2))
grid.draw(v2)
grid.text(bquote(bold(""~"Common")), y=0.93, gp = gpar(cex=1.2))
grid.text(LETTERS[2], y=0.93, x=0.07, gp = gpar(cex=1.6))
popViewport()
dev.off()
