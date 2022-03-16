library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)

# loading distal_only
load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/Distal_Only_Enhancers.Rda")
# loading p5k (promoter contacts over 5kb away)
load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/over5k_contacts.Rda")
# getting expression functions .maxExpressionInsideRegions and .medianExpressionInsideRegions
source("/home/jw18713/Project1/Paper_Source_Scripts/expression_functions_BG3.R")
# expression data for BG3
gene_expression_data <- read.csv("/home/jw18713/Project1/data/gene_expression_data.csv")

distal_only <- distal_only[width(distal_only) >= 1000] # Just checking super enhancers

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
list_distal_only <- split(distal_only, as.factor(distal_only))

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
distal_only_glist <- lapply(list_distal_only, expression_extraction, expression_data = contacted_expressed)
new <- Sys.time()
print(old - new)
save(distal_only_glist, file = "/home/jw18713/Project1/Paper_Rda_Stuff/BG3/distal_only_glist_SE.Rda")

# Now which of these don't have any expresion?
no_expression <- distal_only_glist[sapply(distal_only_glist, length) == 0]

# And which of these are expressed?
distal_only_contacted <- distal_only_glist[sapply(distal_only_glist, length) > 0]

# Hey, this works!
# Now we need median expression data
medianExp <- function(gr){
  return(median(gr$BG3_FPKM, na.rm = T))
}

medians_plot <- sapply(distal_only_contacted, medianExp)

# And I may as well sort the max one as well
maxExp <- function(gr){
  return(max(gr$BG3_FPKM, na.rm = T))
}

max_to_plot <- sapply(distal_only_contacted, maxExp)

no_expression_median <- c(no_expression, distal_only_contacted[medians_plot == 0])
no_expression_max <- c(no_expression, distal_only_contacted[max_to_plot == 0])
expression_median <- medians_plot[medians_plot > 0]
expression_max <- max_to_plot[max_to_plot > 0]

# I can plot this for both methods below
pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure4/BG3/BG3_Distal_Super_Enhancers_Zero_Expression_Barplot_Median.pdf")
barplot(c(length(expression_median), length(no_expression_median)),
      main = "Contacted Regions with Expression",
      xlab = "Expression of Contacted Regions",
      ylab = "Number of Predicted Enhancers",
      col = "lightblue",
      names.arg = c("Expressed", "Not Expressed"),
      ylim = c(0,1500))
dev.off()

pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure4/BG3/BG3_Distal_Super_Enhancers_Zero_Expression_Barplot_Max.pdf")
barplot(c(length(expression_max), length(no_expression_max)),
      main = "Contacted Regions with Expression",
      xlab = "Expression of Contacted Regions",
      ylab = "Number of Predicted Enhancers",
      col = "lightblue",
      names.arg = c("Expressed", "Not Expressed"),
      ylim = c(0,1500))
dev.off()


# Ok, so this is my histogram then
pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure4/BG3/BG3_Distal_Super_Enhancer_Expression_Histogram.pdf")
hist(log2(expression_max+1), main="BG3 Distal Only Super Enhancers\nMax Contacted Expression", xlab="Expression FPKM", col="lightblue", cex.main=2, cex.lab=1.5, axes=F, ylim = c(0,500))
axis(1, at=seq(0,15,by=1), labels=c(0,2^seq(1,15, by = 1)))
axis(2, at=seq(0,500, by = 100), labels = seq(0,500, by = 100))
dev.off()
