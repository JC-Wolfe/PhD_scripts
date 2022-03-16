library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(VennDiagram)

load("/home/jw18713/Archive_Year1/grown_predicted_enhancers_S2.Rda") # pred_enh
gr_contacts <- get(load("/home/jw18713/Archive_Year1/useful_HiC_dm6_contactMap_S2.Rda"))
ann_build <- get(load(file = "/home/jw18713/annotated_GRange_FI_done.Rda"))
gene_expression_data <- read.csv("/home/jw18713/Project1/data/gene_expression_data.csv")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
load("/home/jw18713/Archive_Year1/starr_S2.Rda") # as starr

shared_enhancers <- pred_enh[pred_enh %over% starr]

gr_contacts$contact_chr <- paste0("chr", as.character(gr_contacts$contact_chr))

contact_over5k <- gr_contacts[end(gr_contacts) + 5000 < gr_contacts$contact_start
                              |!as.character(seqnames(gr_contacts)) == gr_contacts$contact_chr]

p5k <- contact_over5k[contact_over5k$promoter_contact == TRUE] # Contact promoter over 5kb from enhancer

contact_under5k <- gr_contacts[end(gr_contacts) + 5000 > gr_contacts$contact_start
                              & as.character(seqnames(gr_contacts)) == gr_contacts$contact_chr]
p_under5k <- contact_under5k[contact_under5k$promoter_contact == TRUE] # Contact promoter under 5kb from enhancer

# Distal stuff
starr_overlaps <- findOverlaps(starr, p5k) # Distal STARR Overlaps
distal_enhancer_overlaps <- findOverlaps(pred_enh, p5k) # Distal Enhancer Overlaps
shared_enhancer_overlaps <- findOverlaps(shared_enhancers, p5k)

starr_prom_contacts <- reduce(starr[queryHits(starr_overlaps)]) # Distal STARR with promoter contacts
distal_enhancers <- reduce(pred_enh[queryHits(distal_enhancer_overlaps)]) # Distal Enhancers with promoter contacts
shared_distal_enhancers <- reduce(shared_enhancers[queryHits(shared_enhancer_overlaps)])

# Proximal stuff
starr_prox_overlaps <- findOverlaps(starr, p_under5k) # Proximal STARR overlaps
proximal_enhancer_overlaps <- findOverlaps(pred_enh, p_under5k) # Proximal Enhancer overlaps
proximal_shared <- findOverlaps(shared_enhancers, p_under5k)

starr_prox_contacts <- reduce(starr[queryHits(starr_prox_overlaps)]) # Proximal STARR with promoter contacts
proximal_enhancers <- reduce(pred_enh[queryHits(proximal_enhancer_overlaps)]) # Proximal enhancers with promoter contacts
shared_proximal_enhancers <- reduce(shared_enhancers[queryHits(proximal_shared)])

# Neither enhancers
starr_neither <- starr[!starr %over% starr_prom_contacts & !starr %over% starr_prox_contacts] # Do not overlap a promoter
enhancers_neither <- pred_enh[!pred_enh %over% distal_enhancers & !pred_enh %over% proximal_enhancers] # Do not overlap a promoter
shared_neither <- shared_enhancers[!shared_enhancers %over% shared_distal_enhancers & !shared_enhancers %over% shared_proximal_enhancers]

# Distal Only
distal_only <- distal_enhancers[!distal_enhancers %over% proximal_enhancers]
shared_distal_only <- shared_distal_enhancers[!shared_distal_enhancers %over% shared_proximal_enhancers]



save(enhancers_neither, file = "/home/jw18713/Project1/Paper_Rda_Stuff/S2/Pred_Enh_No_Promoter_Overlaps.Rda")
save(starr_neither, file = "/home/jw18713/Project1/Paper_Rda_Stuff/S2/Starr_No_Promoter_Overlaps.Rda")
save(proximal_enhancers, file = "/home/jw18713/Project1/Paper_Rda_Stuff/S2/Proximal_Enhancers.Rda")
save(distal_only, file = "/home/jw18713/Project1/Paper_Rda_Stuff/S2/Distal_Only_Enhancers.Rda")
save(distal_enhancers, file = "/home/jw18713/Project1/Paper_Rda_Stuff/S2/All_Distal_Enhancers.Rda")
save(starr, file = "/home/jw18713/Project1/Paper_Rda_Stuff/S2/starr.Rda")
save(p5k, file = "/home/jw18713/Project1/Paper_Rda_Stuff/S2/over5k_contacts.Rda")
save(p_under5k, file = "/home/jw18713/Project1/Paper_Rda_Stuff/S2/under5k_contacts.Rda")
save(shared_distal_enhancers, file = "/home/jw18713/Project1/Paper_Rda_Stuff/S2/shared_distal_enhancers.Rda")
save(shared_proximal_enhancers, file = "/home/jw18713/Project1/Paper_Rda_Stuff/S2/shared_proximal_enhancers.Rda")
save(shared_neither, file = "/home/jw18713/Project1/Paper_Rda_Stuff/S2/shared_neither.Rda")
save(shared_distal_only, file = "/home/jw18713/Project1/Paper_Rda_Stuff/S2/shared_distal_only.Rda")


# Right, so now we have 4 things to work with realisitically
# We have:
# Distal STARR overlaps
# Distal Enhancer overlaps
# Proximal STARR overlaps
# Proximal Enhancer overlaps


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
pdf("/home/jw18713/Project1/Paper_Plots/Figure3/S2/S2_Enhancer_VennDiagram.pdf")
plot.new()
v1 <- generateTripleVenn(distal_enhancers, proximal_enhancers, enhancers_neither, c("Distal\nPromoter", "Proximal\nPromoter", "Neither"),
                           cols=c(cbbPalette[6], cbbPalette[7], cbbPalette[8]),
                           cat.pos = c(-60, 60, 180),
                           cat.dist = c(0.125, 0.125, 0.05),
                           human = F)
title(main = "S2 Novel Enhancer/Promoter Overlaps")
dev.off()


# Distal Enhancers
pdf("/home/jw18713/Project1/Paper_Plots/Figure3/S2/S2_Distal_Only_Enhancer_Widths.pdf")
hist(log2(width(distal_only)), axes=F, main="S2 Novel Distal Only Enhancer Widths", xlab="Width (bp)", col="lightblue", cex.main=2, cex.lab=1.5, xlim = c(3,14), ylim = c(0,1500))
axis(1, at=seq(3,14, by = 1), labels=2^seq(3,14,by=1))
axis(2, at=seq(0,1500, by = 250), labels = seq(0,1500, by = 250))
abline(v = log2(50), col="red", lty=2, lwd=2)
abline(v = log2(1000), col="red", lty=2, lwd=2)
text(4.25, 1425, "Fragments")
text(8.25, 1425, "Enhancers")
text(12, 1425, "Super\nEnhancers")
dev.off()

# Proximal Enhancer
pdf("/home/jw18713/Project1/Paper_Plots/Figure3/S2/S2_Proximal_Enhancer_Widths.pdf")
hist(log2(width(proximal_enhancers)), axes=F, main="S2 Novel Proximal Enhancer Widths", xlab="Width (bp)", col="lightblue", cex.main=2, cex.lab=1.5, xlim = c(3,14), ylim = c(0,1500))
axis(1, at=seq(3,14, by = 1), labels=2^seq(3,14,by=1))
axis(2, at=seq(0,1500, by = 250), labels = seq(0,1500, by = 250))
abline(v = log2(50), col="red", lty=2, lwd=2)
abline(v = log2(1000), col="red", lty=2, lwd=2)
text(4.25, 1425, "Fragments")
text(8.25, 1425, "Enhancers")
text(12, 1425, "Super\nEnhancers")
dev.off()
