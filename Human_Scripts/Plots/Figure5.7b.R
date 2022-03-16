# VennDiagram package
if(!require("VennDiagram", character.only = TRUE)){
  install.packages("VennDiagram")
}
library(VennDiagram)

# install package gridExtra if not already installed
if(!require("gridExtra", character.only = TRUE)){
  install.packages("gridExtra")
}
library(gridExtra)

# install package grid if not already installed
if(!require("grid", character.only = TRUE)){
  install.packages("grid")
}
library(grid)

# install package gridBase if not already installed
if(!require("gridBase", character.only = TRUE)){
  install.packages("gridBase")
}
library(gridBase)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

generateTripleVenn <- function(set1, set2, set3, categories, cols=c("#0072B2", "#D55E00", "#CC79A7"), cat.pos = c(-20, 0, 20), cat.dist = c(0.07, 0.07, 0.07), human = TRUE){

  v <- draw.triple.venn(area1 = length(set1), area2 = length(set2), area3 = length(set3),
                          n12 = sum(set1%over%set2),
                          n23 = sum(set2%over%set3),
                          n13 = sum(set1%over%set3),
                          n123 = sum(set1%over%set2 & set1%over%set3),
                          category = categories, col = "transparent", fill = cols,
                          alpha = 0.6, label.col = rep("black", 7), cex = 1.2,
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

# Loading Enhancer Atlas
EA <- import("~/Human_Datasets/H9/E_Atlas/H9.bed")
E_Atlas <- reduce(EA, min.gapwidth = 1000)

# Loading H9
load("~/Human_Datasets/H9/Processed_datasets/H9_XAI.Rda")
H9 <- H9_XAI[H9_XAI$Conf_Perc_1 >= 0.8]
H9 <- reduce(H9, min.gapwidth = 1000)

STARR <- reduce(H9_XAI[H9_XAI$STARR_seq_binary == 1], min.gapwidth = 1000)

labels <- c(paste0("Enhancer\nAtlas 2.0\n(", length(E_Atlas), ")"),
  paste0("H9\n(", length(H9), ")"),
  paste0("STARR-seq\n(", length(STARR), ")"))

v1 <- generateTripleVenn(E_Atlas, H9, STARR, c(labels),
                           cols=c(cbbPalette[6], cbbPalette[7], cbbPalette[8]),
                           cat.pos = c(-60, 60, 180),
                           cat.dist = c(0.125, 0.125, 0.125),
                           human = F)

dev.off()

pdf("/home/jw18713/Human_Datasets/Plots/Supp/EA_H9_STARR-seq_VD.pdf", width=9, height=10,pointsize = 14);
par(mar=c(0, 0, 0, 0)+0.1)
pushViewport(plotViewport(layout=grid.layout(1, 1),gp=gpar(cex=1.0)))
pushViewport(plotViewport(layout.pos.row=1, layout.pos.col=1))
grid.draw(v1)
grid.text(bquote(bold("Predicted Enhancers")), y=0.93, gp = gpar(cex=1.2))
popViewport()
dev.off()
