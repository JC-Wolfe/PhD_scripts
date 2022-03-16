library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(VennDiagram)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/Pred_Enh_No_Promoter_Overlaps.Rda") #enhancers_neither
load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/Proximal_Enhancers.Rda") # proximal_enhancers
load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/Distal_Only_Enhancers.Rda") # distal_only
load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/All_Distal_Enhancers.Rda") # distal_enhancers

distal_only_SE <- distal_only[width(distal_only) >= 1000]
proximal_enhancers_SE <- proximal_enhancers[width(proximal_enhancers) >= 1000]
enhancers_neither_SE <- enhancers_neither[width(enhancers_neither) >= 1000]
distal_enhancers_SE <- distal_enhancers[width(distal_enhancers) >= 1000]

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

# VENN DIAGRAMS
# Plotting the predicted enhancer locations
pdf("/home/jw18713/Project1/Paper_Plots/Figure4/S2/S2_Super_Enhancer_VennDiagram.pdf")
plot.new()
v1 <- generateTripleVenn(distal_enhancers_SE, proximal_enhancers_SE, enhancers_neither_SE, c("Distal\nPromoter", "Proximal\nPromoter", "Neither"),
                           cols=c(cbbPalette[6], cbbPalette[7], cbbPalette[8]),
                           cat.pos = c(-60, 60, 180),
                           cat.dist = c(0.125, 0.125, 0.05),
                           human = F)
title(main = "S2 Novel Super Enhancer/Promoter Overlaps")
dev.off()


# Distal Enhancers
pdf("/home/jw18713/Project1/Paper_Plots/Figure4/S2/S2_Distal_Only_Super_Enhancer_Widths.pdf")
hist(width(distal_only_SE), axes=F, main="S2 Novel Distal Only\nSuper Enhancer Widths", xlab="Width (bp)", col="lightblue", cex.main=2, cex.lab=1.5, xlim = c(1000,15000), ylim = c(0,1000))
axis(1, at=seq(1000,11000, by = 1000), labels=seq(1000,11000,by=1000))
axis(2, at=seq(0,1000, by = 100), labels = seq(0,1000, by = 100))
dev.off()

# Proximal Enhancer
pdf("/home/jw18713/Project1/Paper_Plots/Figure4/S2/S2_Proximal_Super_Enhancer_Widths.pdf")
hist(width(proximal_enhancers_SE), axes=F, main="S2 Novel Proximal\nSuper Enhancer Widths", xlab="Width (bp)", col="lightblue", cex.main=2, cex.lab=1.5, xlim = c(1000,15000), ylim = c(0,1000))
axis(1, at=seq(1000,15000, by = 1000), labels=seq(1000,15000,by=1000))
axis(2, at=seq(0,1000, by = 100), labels = seq(0,1000, by = 100))
dev.off()
