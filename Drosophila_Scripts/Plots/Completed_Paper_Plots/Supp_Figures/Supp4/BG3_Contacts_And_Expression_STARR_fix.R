library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(VennDiagram)

gr_contacts <- get(load("~/Archive_Year1/useful_HiC_dm6_contactMap.Rda"))
ann_build <- get(load(file = "/home/jw18713/Archive_Year1/annotated_GRange_FI_done.Rda"))
gene_expression_data <- read.csv("/home/jw18713/Project1/data/gene_expression_data.csv")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
starr <- get(load("~/Project1/Paper_Rda_Stuff/BG3/starr.Rda")) # as starr

gr_contacts$contact_chr <- paste0("chr", as.character(gr_contacts$contact_chr))

contact_over5k <- gr_contacts[end(gr_contacts) + 5000 < gr_contacts$contact_start
                              |!as.character(seqnames(gr_contacts)) == gr_contacts$contact_chr]
p5k <- contact_over5k[contact_over5k$promoter_contact == TRUE] # Contact promoter over 5kb from enhancer

contact_under5k <- gr_contacts[end(gr_contacts) + 5000 > gr_contacts$contact_start
                              & as.character(seqnames(gr_contacts)) == gr_contacts$contact_chr]
p_under5k <- contact_under5k[contact_under5k$promoter_contact == TRUE] # Contact promoter under 5kb from enhancer

# Distal stuff
starr_overlaps <- findOverlaps(starr, p5k) # Distal STARR Overlaps
distal_starr <- reduce(starr[queryHits(starr_overlaps)]) # Distal STARR with promoter contacts

# Proximal stuff
starr_prox_overlaps <- findOverlaps(starr, p_under5k) # Proximal STARR overlaps
proximal_starr <- reduce(starr[queryHits(starr_prox_overlaps)]) # Proximal STARR with promoter contacts

# Neither enhancers
starr_neither <- starr[!starr %over% distal_starr & !starr %over% proximal_starr] # Do not overlap a promoter

# Distal Only
starr_distal_only <-distal_starr[!distal_starr %over% proximal_starr]



save(starr_neither, file = "/home/jw18713/Project1/Paper_Rda_Stuff/BG3/starr_neither.Rda")
save(distal_starr, file = "/home/jw18713/Project1/Paper_Rda_Stuff/BG3/distal_starr.Rda")
save(proximal_starr, file = "/home/jw18713/Project1/Paper_Rda_Stuff/BG3/proximal_starr.Rda")
save(starr_distal_only, file = "/home/jw18713/Project1/Paper_Rda_Stuff/BG3/starr_distal_only.Rda")
