library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(gplots)
library(lattice)

S2_enhancers <- get(load("Archive_Year1/grown_predicted_enhancers_S2.Rda"))
S2_starr <- get(load("Archive_Year1/starr_S2.Rda"))

S2_fuzzy <- S2_enhancers[!S2_enhancers %over% S2_starr]

S2_novel <- S2_fuzzy

load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/All_Distal_Enhancers.Rda") #distal_enhancers
load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/Distal_Only_Enhancers.Rda") #distal_only
# loading proximal_enhancers
load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/Proximal_Enhancers.Rda")

S2_novel_proximal <- S2_novel[S2_novel %over% proximal_enhancers]
S2_novel_distal <- S2_novel[S2_novel %over% distal_only]
S2_novel_neither <- S2_novel[!S2_novel %over% proximal_enhancers & !S2_novel %over% distal_only]
S2_novel_distal_all <- S2_novel[S2_novel %over% distal_enhancers]

save(S2_novel, file = "Archive_Year1/S2_novel.Rda")
save(S2_novel_neither, file = "Archive_Year1/S2_novel_neither.Rda")
save(S2_novel_proximal, file = "Archive_Year1/S2_novel_proximal.Rda")
save(S2_novel_distal, file = "Archive_Year1/S2_novel_distal.Rda")
save(S2_novel_distal_all, file = "Archive_Year1/S2_novel_distal_all.Rda")
