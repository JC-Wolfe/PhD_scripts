
library(GenomicRanges)
library(rtracklayer)

# loading distal_only
load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/background_distal_only.Rda")
# loading p5k (promoter contacts over 5kb away)
load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/over5k_contacts.Rda")
# getting expression functions .maxExpressionInsideRegions and .medianExpressionInsideRegions
source("/home/jw18713/Project1/Paper_Source_Scripts/expression_functions_BG3.R")
# expression data for BG3
gene_expression_data <- read.csv("/home/jw18713/Project1/data/gene_expression_data.csv")

# I need to do an expression bar chart and histogram, for now I'll do this for max and median
# These are the median plots

expression <- GRanges(seqnames = gene_expression_data$Chromosome.arm, ranges=IRanges(gene_expression_data$Start, gene_expression_data$End))
expression$Gene_short_name <- gene_expression_data$Gene.short.name
expression$BG3_FPKM <- gene_expression_data$BG3_FPKM

# So lets see where we have promoter contacts then
exp_test <- GRanges(seqnames = p5k$contact_chr, ranges = IRanges(p5k$contact_start, p5k$contact_end))
exp_over <- findOverlaps(exp_test, expression)
contacted_expressed <- p5k[queryHits(exp_over)]
contacted_expressed$Gene_short_name <- expression$Gene_short_name[subjectHits(exp_over)]
contacted_expressed$BG3_FPKM <- expression$BG3_FPKM[subjectHits(exp_over)]

# This is a GRange list that will contain all of the target enhancers position by position
list_distal_only <- split(background_distal_only, as.factor(background_distal_only))

expression_extraction <- function(gr, expression_data){
  ovr <- findOverlaps(gr, expression_data)
  gr <- gr[queryHits(ovr)]
  gr$contact_chr <- expression_data[subjectHits(ovr)]$contact_chr
  gr$contact_start <- expression_data[subjectHits(ovr)]$contact_start
  gr$contact_end <- expression_data[subjectHits(ovr)]$contact_end
  gr$Gene_short_name <- expression_data[subjectHits(ovr)]$Gene_short_name
  gr$BG3_FPKM <- expression_data[subjectHits(ovr)]$BG3_FPKM
  return(gr)
}

old <- Sys.time()
BG3_background_distal_only_glist <- lapply(list_distal_only, expression_extraction, expression_data = contacted_expressed)
new <- Sys.time()
print(old - new)
save(BG3_background_distal_only_glist, file = "/home/jw18713/Project1/Paper_Rda_Stuff/BG3/BG3_background_distal_only_glist.Rda")

rm(list=ls())
# Everything from here down is for proximal enhancers

# loading proximal_enhancers
load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/proximal_background_contacts.Rda")
# loading p_under5k (promoter contacts under 5kb away)
load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/under5k_contacts.Rda")
# getting expression functions .maxExpressionInsideRegions and .medianExpressionInsideRegions
source("/home/jw18713/Project1/Paper_Source_Scripts/expression_functions_BG3.R")
# expression data for BG3
gene_expression_data <- read.csv("/home/jw18713/Project1/data/gene_expression_data.csv")

# I need to do an expression bar chart and histogram, for now I'll do this for max and median
# These are the median plots

expression <- GRanges(seqnames = gene_expression_data$Chromosome.arm, ranges=IRanges(gene_expression_data$Start, gene_expression_data$End))
expression$Gene_short_name <- gene_expression_data$Gene.short.name
expression$BG3_FPKM <- gene_expression_data$BG3_FPKM

# So lets see where we have promoter contacts then
exp_test <- GRanges(seqnames = p_under5k$contact_chr, ranges = IRanges(p_under5k$contact_start, p_under5k$contact_end))
exp_over <- findOverlaps(exp_test, expression)
contacted_expressed <- p_under5k[queryHits(exp_over)]
contacted_expressed$Gene_short_name <- expression$Gene_short_name[subjectHits(exp_over)]
contacted_expressed$BG3_FPKM <- expression$BG3_FPKM[subjectHits(exp_over)]

# This is a GRange list that will contain all of the target enhancers position by position
list_proximal_enhancers <- split(proximal_background_contacts, as.factor(proximal_background_contacts))

expression_extraction <- function(gr, expression_data){
  ovr <- findOverlaps(gr, expression_data)
  gr <- gr[queryHits(ovr)]
  gr$contact_chr <- expression_data[subjectHits(ovr)]$contact_chr
  gr$contact_start <- expression_data[subjectHits(ovr)]$contact_start
  gr$contact_end <- expression_data[subjectHits(ovr)]$contact_end
  gr$Gene_short_name <- expression_data[subjectHits(ovr)]$Gene_short_name
  gr$BG3_FPKM <- expression_data[subjectHits(ovr)]$BG3_FPKM
  return(gr)
}

old <- Sys.time()
BG3_proximal_background_glist <- lapply(list_proximal_enhancers, expression_extraction, expression_data = contacted_expressed)
new <- Sys.time()
print(old - new)
save(BG3_proximal_background_glist, file = "/home/jw18713/Project1/Paper_Rda_Stuff/BG3/BG3_proximal_background_glist.Rda")

rm(list=ls())

# Everything below here is for S2

# loading distal_only
load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/background_distal_only.Rda")
# loading p5k (promoter contacts over 5kb away)
load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/over5k_contacts.Rda")
# getting expression functions .maxExpressionInsideRegions and .medianExpressionInsideRegions
source("/home/jw18713/Project1/Paper_Source_Scripts/expression_functions_S2.R")
# expression data for S2
gene_expression_data <- read.csv("/home/jw18713/Project1/data/gene_expression_data.csv")


# I need to do an expression bar chart and histogram, for now I'll do this for max and median
# These are the median plots

expression <- GRanges(seqnames = gene_expression_data$Chromosome.arm, ranges=IRanges(gene_expression_data$Start, gene_expression_data$End))
expression$Gene_short_name <- gene_expression_data$Gene.short.name
expression$S2_FPKM <- gene_expression_data$S2DRSC_FPKM

# So lets see where we have promoter contacts then
exp_test <- GRanges(seqnames = p5k$contact_chr, ranges = IRanges(p5k$contact_start, p5k$contact_end))
exp_over <- findOverlaps(exp_test, expression)
contacted_expressed <- p5k[queryHits(exp_over)]
contacted_expressed$Gene_short_name <- expression$Gene_short_name[subjectHits(exp_over)]
contacted_expressed$S2_FPKM <- expression$S2_FPKM[subjectHits(exp_over)]

# This is a GRange list that will contain all of the target enhancers position by position
list_distal_only <- split(background_distal_only, as.factor(background_distal_only))

expression_extraction <- function(gr, expression_data){
  ovr <- findOverlaps(gr, expression_data)
  gr <- gr[queryHits(ovr)]
  gr$contact_chr <- expression_data[subjectHits(ovr)]$contact_chr
  gr$contact_start <- expression_data[subjectHits(ovr)]$contact_start
  gr$contact_end <- expression_data[subjectHits(ovr)]$contact_end
  gr$Gene_short_name <- expression_data[subjectHits(ovr)]$Gene_short_name
  gr$S2_FPKM <- expression_data[subjectHits(ovr)]$S2_FPKM
  return(gr)
}

old <- Sys.time()
S2_background_distal_only_glist <- lapply(list_distal_only, expression_extraction, expression_data = contacted_expressed)
new <- Sys.time()
print(old - new)
save(S2_background_distal_only_glist, file = "/home/jw18713/Project1/Paper_Rda_Stuff/S2/S2_background_distal_only_glist.Rda")

rm(list=ls())

# Everything below here is S2 proximal


# loading proximal_enhancers
load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/proximal_background_contacts.Rda")
# loading p_under5k (promoter contacts under 5kb away)
load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/under5k_contacts.Rda")
# getting expression functions .maxExpressionInsideRegions and .medianExpressionInsideRegions
source("/home/jw18713/Project1/Paper_Source_Scripts/expression_functions_S2.R")
# expression data for S2
gene_expression_data <- read.csv("/home/jw18713/Project1/data/gene_expression_data.csv")

# I need to do an expression bar chart and histogram, for now I'll do this for max and median
# These are the median plots

expression <- GRanges(seqnames = gene_expression_data$Chromosome.arm, ranges=IRanges(gene_expression_data$Start, gene_expression_data$End))
expression$Gene_short_name <- gene_expression_data$Gene.short.name
expression$S2_FPKM <- gene_expression_data$S2DRSC_FPKM

# So lets see where we have promoter contacts then
exp_test <- GRanges(seqnames = p_under5k$contact_chr, ranges = IRanges(p_under5k$contact_start, p_under5k$contact_end))
exp_over <- findOverlaps(exp_test, expression)
contacted_expressed <- p_under5k[queryHits(exp_over)]
contacted_expressed$Gene_short_name <- expression$Gene_short_name[subjectHits(exp_over)]
contacted_expressed$S2_FPKM <- expression$S2_FPKM[subjectHits(exp_over)]

# This is a GRange list that will contain all of the target enhancers position by position
list_proximal_enhancers <- split(proximal_background_contacts, as.factor(proximal_background_contacts))

expression_extraction <- function(gr, expression_data){
  ovr <- findOverlaps(gr, expression_data)
  gr <- gr[queryHits(ovr)]
  gr$contact_chr <- expression_data[subjectHits(ovr)]$contact_chr
  gr$contact_start <- expression_data[subjectHits(ovr)]$contact_start
  gr$contact_end <- expression_data[subjectHits(ovr)]$contact_end
  gr$Gene_short_name <- expression_data[subjectHits(ovr)]$Gene_short_name
  gr$S2_FPKM <- expression_data[subjectHits(ovr)]$S2_FPKM
  return(gr)
}

old <- Sys.time()
S2_proximal_background_glist <- lapply(list_proximal_enhancers, expression_extraction, expression_data = contacted_expressed)
new <- Sys.time()
print(old - new)
save(S2_proximal_background_glist, file = "/home/jw18713/Project1/Paper_Rda_Stuff/S2/S2_proximal_background_glist.Rda")

rm(list=ls())
