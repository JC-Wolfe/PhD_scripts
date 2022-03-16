library(BSgenome.Dmelanogaster.UCSC.dm6)
library(GenomicRanges)
library(rtracklayer)
library(gplots)
library(lattice)

BG3_enhancers <- get(load("Archive_Year1/grown_predicted_enhancers_BG3.Rda"))
BG3_starr <- get(load("Archive_Year1/starr.Rda"))

BG3_fuzzy <- BG3_enhancers[!BG3_enhancers %over% BG3_starr]

BG3_novel <- BG3_fuzzy

load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/Distal_Only_Enhancers.Rda") # distal_only
load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/All_Distal_Enhancers.Rda") # distal_enhancers
# loading proximal_enhancers
load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/Proximal_Enhancers.Rda")

BG3_novel_proximal <- BG3_novel[BG3_novel %over% proximal_enhancers]
BG3_novel_distal <- BG3_novel[BG3_novel %over% distal_only]
BG3_novel_neither <- BG3_novel[!BG3_novel %over% proximal_enhancers & !BG3_novel %over% distal_only]
BG3_novel_distal_all <- BG3_novel[BG3_novel %over% distal_enhancers]

save(BG3_novel, file = "Archive_Year1/BG3_novel.Rda")
save(BG3_novel_neither, file = "Archive_Year1/BG3_novel_neither.Rda")
save(BG3_novel_proximal, file = "Archive_Year1/BG3_novel_proximal.Rda")
save(BG3_novel_distal, file = "Archive_Year1/BG3_novel_distal.Rda")
save(BG3_novel_distal_all, file = "Archive_Year1/BG3_novel_distal_all.Rda")
