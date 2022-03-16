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


BG3_bg_distal <- get(load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/distal_background_contacts.Rda"))
BG3_bg_proximal <- get(load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/proximal_background_contacts.Rda"))
BG3_bg_neither <- get(load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/background_neither.Rda"))

v1 <- generateTripleVenn(BG3_bg_distal, BG3_bg_proximal, BG3_bg_neither, c("Distal\nPromoter", "Proximal\nPromoter", "Neither"),
                           cols=c(cbbPalette[6], cbbPalette[7], cbbPalette[8]),
                           cat.pos = c(-60, 60, 180),
                           cat.dist = c(0.125, 0.125, 0.05),
                           human = F)


S2_bg_distal <- get(load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/distal_background_contacts.Rda"))
S2_bg_proximal <- get(load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/proximal_background_contacts.Rda"))
S2_bg_neither <- get(load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/background_neither.Rda"))




v2 <- generateTripleVenn(S2_bg_distal, S2_bg_proximal, S2_bg_neither, c("Distal\nPromoter", "Proximal\nPromoter", "Neither"),
                           cols=c(cbbPalette[6], cbbPalette[7], cbbPalette[8]),
                           cat.pos = c(-60, 60, 180),
                           cat.dist = c(0.125, 0.125, 0.05),
                           human = F)

pdf("~/Project1/Paper_Plots/Figure3/Background_contact_venn.pdf", width=8, height=5,pointsize = 10);
par(mar=c(0, 0, 0, 0)+0.1)
pushViewport(plotViewport(layout=grid.layout(1, 2),gp=gpar(cex=1.0)))
pushViewport(plotViewport(layout.pos.row=1, layout.pos.col=1))
grid.draw(v1)
grid.text(bquote(bold(""~"BG3")), y=0.93, gp = gpar(cex=1.2))
grid.text(LETTERS[1], y=0.93, x=0.07, gp = gpar(cex=1.6))
popViewport()
pushViewport(plotViewport(layout.pos.row=1, layout.pos.col=2))
grid.draw(v2)
grid.text(bquote(bold(""~"S2")), y=0.93, gp = gpar(cex=1.2))
grid.text(LETTERS[2], y=0.93, x=0.07, gp = gpar(cex=1.6))
popViewport()
dev.off()
