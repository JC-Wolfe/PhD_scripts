library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(VennDiagram)

load("grown_predicted_enhancers_BG3.Rda") # pred_enh
gr_contacts <- get(load("useful_HiC_dm6_contactMap.Rda"))
ann_build <- get(load(file = "/home/jw18713/annotated_GRange_FI_done.Rda"))
gene_expression_data <- read.csv("/home/jw18713/Project1/data/gene_expression_data.csv")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
starr <- load("starr.Rda")


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

starr_prom_contacts <- reduce(starr[queryHits(starr_overlaps)]) # Distal STARR with promoter contacts
distal_enhancers <- reduce(pred_enh[queryHits(distal_enhancer_overlaps)]) # Distal Enhancers with promoter contacts

# Proximal stuff
starr_prox_overlaps <- findOverlaps(starr, p_under5k) # Proximal STARR overlaps
proximal_enhancer_overlaps <- findOverlaps(pred_enh, p_under5k) # Proximal Enhancer overlaps

starr_prox_contacts <- reduce(starr[queryHits(starr_prox_overlaps)]) # Proximal STARR with promoter contacts
proximal_enhancers <- reduce(pred_enh[queryHits(proximal_enhancer_overlaps)]) # Proximal enhancers with promoter contacts

# Neither enhancers
starr_neither <- starr[!starr %over% starr_prom_contacts & !starr %over% starr_prox_contacts]
enhancers_neither <- pred_enh[!pred_enh %over% distal_enhancers & !pred_enh %over% proximal_enhancers]

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


# Plotting the predicted enhancer locations
pdf("/home/jw18713/Project1/Paper_Plots/BG3_Distal_Enhancer_VennDiagram.pdf")
v1 <- generateTripleVenn(distal_enhancers, proximal_enhancers, enhancers_neither, c("Distal\nEnhancers", "Proximal\nEnhancers", "Neither"),
                           cols=c(cbbPalette[6], cbbPalette[7], cbbPalette[8]),
                           cat.pos = c(-60, 60, 180),
                           cat.dist = c(0.125, 0.125, 0.05),
                           human = F)

dev.off()


# Plotting the STARR enhancer locations
pdf("/home/jw18713/Project1/Paper_Plots/BG3_Distal_STARR_VennDiagram.pdf")
v1 <- generateTripleVenn(starr_prom_contacts, starr_prox_contacts, starr_neither, c("Distal\nSTARR-seq", "Proximal\nSTARR-seq", "Neither"),
                           cols=c(cbbPalette[6], cbbPalette[7], cbbPalette[8]),
                           cat.pos = c(-60, 60, 180),
                           cat.dist = c(0.125, 0.125, 0.05),
                           human = F)

dev.off()
