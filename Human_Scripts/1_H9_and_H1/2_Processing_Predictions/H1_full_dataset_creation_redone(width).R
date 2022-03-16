library(GenomicRanges)
library(rtracklayer)
library(genomation)
library(DMRcaller)

load("~/Human_Datasets/H1/Processed_datasets/H1_XAI.Rda")
predicted <- H1_XAI[H1_XAI$Conf_Perc_1 >= 0.8]
predicted <- reduce(predicted, min.gapwidth = 1000)
predicted <- predicted[width(predicted) >= 200]

load("~/Human_Datasets/H1/DMRcaller_norm_tiled.Rda")

# Chain for lift overs
chain <- import.chain("/home/jw18713/Human_Datasets/hg19ToHg38.over.chain")

# Removing the columns we don't use in the model
H9_full$ATAC_signalValues <- NULL

# Fixing the names by removing the underscores and numbers
fixed_names <- unlist(strsplit(names(mcols(H1_tiled)[1:27]), "_"))[c(T,F)]
content <- mcols(H1_tiled)
names(content)[1:27] <- fixed_names
mcols(H1_tiled) <- content

# Annotating predicted regions as predicted regions
H1_tiled$P_enh <- 0
predicted_over <- findOverlaps(H1_tiled, predicted)
H1_tiled$P_enh[queryHits(predicted_over)] <- 1

# Loading in H9 data for classification
load("~/Human_Datasets/H9/Processed_datasets/H9_full(widths).Rda")

# Classifying enhancers as Neither, Common, Putative, and STARR-seq only
# Neither unless otherwise specified
H1_tiled$classification <- "Neither"

# Putative
H1_tiled$classification[H1_tiled$P_enh == 1
  & H9_full$P_enh == 0] <- "H1 Specific"

# Common
H1_tiled$classification[H1_tiled$P_enh == 1
  & H9_full$P_enh == 1] <- "Common"

# STARR-seq Only
H1_tiled$classification[H1_tiled$P_enh == 0
  & H9_full$P_enh == 1] <- "H9 Specific"

H1_tiled$classification <- as.factor(H1_tiled$classification)

H1_full <- H1_tiled

# Getting annotations
gtf <- import("~/Human_Datasets/Homo_sapiens.GRCh37.87.chr.gtf")
seqlevelsStyle(gtf) <- "UCSC"
gtf <- unlist(liftOver(gtf, chain))

#priority list
#promoter
#first intron
#other intron
#exon
#5' utr
#3' utr
#intergenic

# Creating and filling the annotation meta column
H1_full$annotations <- "intergenic"

transcript <- gtf[gtf$type == "transcript"]
exons <- gtf[gtf$type == "exon"]
fiveprime <- gtf[gtf$type == "five_prime_utr"]
threeprime <- gtf[gtf$type == "three_prime_utr"]
sense_promoter <- gtf[gtf$type == "gene" & strand(gtf) == "+"]
start(sense_promoter) <- start(sense_promoter) - 250
end(sense_promoter) <- start(sense_promoter) + 249
antisense_promoter <- gtf[gtf$type == "gene" & strand(gtf) == "-"]
end(antisense_promoter) <- end(antisense_promoter) + 250
start(antisense_promoter) <- end(antisense_promoter) - 249
promoters <- reduce(c(sense_promoter, antisense_promoter))

overlaps1 <- findOverlaps(H1_full, transcript) #Finding where the transcript annotations are
H1_full$annotations[queryHits(overlaps1)] <- "transcript"
overlaps2 <- findOverlaps(H1_full, exons)
H1_full$annotations[queryHits(overlaps2)] <- "exon"
overlaps3 <- findOverlaps(H1_full, fiveprime)
H1_full$annotations[H1_full$annotations == "transcript"] <- "intron"
H1_full$annotations[queryHits(overlaps3)] <- "five_prime_utr"
overlaps4 <- findOverlaps(H1_full, threeprime)
H1_full$annotations[queryHits(overlaps4)] <- "three_prime_utr"
overlaps5 <- findOverlaps(H1_full, promoters)
H1_full$annotations[queryHits(overlaps5)] <- "promoter"

sense_genes <- gtf[gtf$type == "gene" & strand(gtf) == "+"]
antisense_genes <- gtf[gtf$type == "gene" & strand(gtf) == "-"]

H1_full$sense_gene <- "NA"
sense_over <- findOverlaps(H1_full, sense_genes)
H1_full$sense_gene[queryHits(sense_over)] <- sense_genes$gene_name[subjectHits(sense_over)]

H1_full$antisense_gene <- "NA"
antisense_over <- findOverlaps(H1_full, antisense_genes)
H1_full$antisense_gene[queryHits(antisense_over)] <- antisense_genes$gene_name[subjectHits(antisense_over)]

save(H1_full, file = "~/Human_Datasets/H1/Processed_datasets/H1_full(widths).Rda")
