library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(VennDiagram)

load("/home/jw18713/Archive_Year1/background_sample.Rda") # background_sample
gr_contacts <- get(load("/home/jw18713/Archive_Year1/useful_HiC_dm6_contactMap_S2.Rda"))
ann_build <- get(load(file = "/home/jw18713/Archive_Year1/annotated_GRange_FI_done.Rda"))

gr_contacts$contact_chr <- paste0("chr", as.character(gr_contacts$contact_chr))

contact_over5k <- gr_contacts[end(gr_contacts) + 5000 < gr_contacts$contact_start
                              |!as.character(seqnames(gr_contacts)) == gr_contacts$contact_chr]

p5k <- contact_over5k[contact_over5k$promoter_contact == TRUE] # Contact promoter over 5kb from enhancer

contact_under5k <- gr_contacts[end(gr_contacts) + 5000 > gr_contacts$contact_start
                              & as.character(seqnames(gr_contacts)) == gr_contacts$contact_chr]
p_under5k <- contact_under5k[contact_under5k$promoter_contact == TRUE] # Contact promoter under 5kb from enhancer

# Distal stuff
distal_background_overlaps <- findOverlaps(background_sample, p5k) # Distal background overlaps
distal_background_contacts <- reduce(background_sample[queryHits(distal_background_overlaps)]) # Distal background with promoter contacts

# Proximal stuff
proximal_background_overlaps <- findOverlaps(background_sample, p_under5k) # Proximal STARR overlaps
proximal_background_contacts <- reduce(background_sample[queryHits(proximal_background_overlaps)]) # Proximal STARR with promoter contacts

# Neither enhancers
background_neither <- background_sample[!background_sample %over% distal_background_contacts & !background_sample %over% proximal_background_contacts] # Do not overlap a promoter

# Distal Only
background_distal_only <- distal_background_contacts[!distal_background_contacts %over% proximal_background_contacts]

save(distal_background_contacts, file = "/home/jw18713/Project1/Paper_Rda_Stuff/S2/distal_background_contacts.Rda")
save(proximal_background_contacts, file = "/home/jw18713/Project1/Paper_Rda_Stuff/S2/proximal_background_contacts.Rda")
save(background_neither, file = "/home/jw18713/Project1/Paper_Rda_Stuff/S2/background_neither.Rda")
save(background_distal_only, file = "/home/jw18713/Project1/Paper_Rda_Stuff/S2/background_distal_only.Rda")
