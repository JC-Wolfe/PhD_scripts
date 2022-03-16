library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
# Getting the hg38 genome
genome <- getSeq(BSgenome.Hsapiens.UCSC.hg38)

# Getting my own functions
source("~/Genome_construction_functions.R")

# Loading in contacts and H9 full genome
setwd("~/Human_Datasets/H9/HiC/Jareth_Liv_HiC2/compartments/jarethTxtFiles")
load("~/Human_Datasets/H9/Processed_datasets/H9_full(widths).Rda")

gr_contacts <- GRanges()

cbin_in <- function(x, chr){
  all <- read.csv(x, sep = "\t")
  flt <- all[!all[,3]=="NaN",] # removing NaNs
  gr <- GRanges(seqnames = chr,
    ranges = IRanges(start = flt[,1], end = flt[,1] + 10000))

  gr$contact_chr <- chr
  gr$contact_start <- flt[,2]
  gr$contact_end <- flt[,2] + 10000

  return(gr)
}

gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr1.txt", chr = "chr1"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr2.txt", chr = "chr2"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr3.txt", chr = "chr3"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr4.txt", chr = "chr4"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr5.txt", chr = "chr5"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr6.txt", chr = "chr6"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr7.txt", chr = "chr7"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr8.txt", chr = "chr8"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr9.txt", chr = "chr9"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr10.txt", chr = "chr10"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr11.txt", chr = "chr11"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr12.txt", chr = "chr12"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr13.txt", chr = "chr13"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr14.txt", chr = "chr14"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr15.txt", chr = "chr15"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr16.txt", chr = "chr16"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr17.txt", chr = "chr17"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr18.txt", chr = "chr18"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr19.txt", chr = "chr19"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr20.txt", chr = "chr20"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr21.txt", chr = "chr21"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chr22.txt", chr = "chr22"))
gr_contacts <- c(gr_contacts, cbin_in("GSM3262957_chrX.txt", chr = "chrX"))


contact_check <- GRanges(seqnames = gr_contacts$contact_chr,
                          ranges = IRanges(gr_contacts$contact_start,
                          gr_contacts$contact_end))
seqlevelsStyle(contact_check) <- "UCSC"

promoters <- reduce(H9_full[H9_full$annotations == "promoter"])
prm <- rep(F, length(contact_check))

overlaps <- findOverlaps(promoters, contact_check)
prm[subjectHits(overlaps)] <- T

gr_contacts$promoter_contact <- prm

save(gr_contacts,
  file = "~/Pre_viva/Hi-C/gr_contacts.Rda")

#------------------------------Contact Mapping---------------------------------#
#------------------------------------------------------------------------------#

contact_mapping <- function(x, contacts, promoters, distance = 5000, threshold = 0){

  over5k <- contacts[(end(contacts) + distance < contacts$contact_start
                    |!as.character(seqnames(contacts)) == contacts$contact_chr)
                    #& contacts$enrichment > threshold
                    & contacts$promoter_contact == TRUE]

  nearest_promoter <- nearest(x, promoters)

  prox <- x[distance(x,promoters[nearest_promoter]) <= distance]

  # Distal overlaps
  dist <- x[x %over% over5k]
  # Proximal
  # Distal and Proximal
  d_and_p <- dist[dist %over% prox]
  # Distal only
  d_only <- dist[!dist %over% prox]
  # Proximal only
  p_only <- prox[!prox %over% dist]
  # Neither
  neither <- x[!x %over% dist & !x %over% prox]

  # GRange List Output
  out <- GRangesList(d_and_p, d_only, p_only, neither)
  names(out) <- c("Distal_&_Proximal_Overlaps",
                  "Distal_Only_Overlaps",
                  "Proximal_Only_Overlaps",
                  "No_Promoter_Overlaps")
  return(out)
}

putative <- reduce(H9_full[H9_full$classification == "Putative"])
common <- reduce(H9_full[H9_full$classification == "Common"])
starr_only <- reduce(H9_full[H9_full$classification == "STARR"])
neither <- reduce(H9_full[H9_full$classification == "Neither"])
background <- empty_tiled_genome(Hsapiens, binsize = 2000, chr_range = seq(1,23))

H9_putative <- contact_mapping(putative, gr_contacts, promoters)
H9_common <- contact_mapping(common, gr_contacts, promoters)
H9_starr_only <- contact_mapping(starr_only, gr_contacts, promoters)
H9_neither <- contact_mapping(neither, gr_contacts, promoters)
H9_background <- contact_mapping(background, gr_contacts, promoters)

setwd("~/Pre_viva/Hi-C")
save(H9_putative, file = "contact_Rda/H9_putative.Rda")
save(H9_common, file = "contact_Rda/H9_common.Rda")
save(H9_starr_only, file = "contact_Rda/H9_starr_ony.Rda")
save(H9_neither, file = "contact_Rda/H9_neither.Rda")
save(H9_background, file = "contact_Rda/H9_background.Rda")

#------------------------------Contact Objects---------------------------------#
#------------------------------------------------------------------------------#

distal_contact_creation <- function(contacts, promoters, distance = 5000,
  threshold = 0){
    distal <-contacts[(end(contacts) + distance < contacts$contact_start
            |!as.character(seqnames(contacts)) == contacts$contact_chr)
            #& contacts$enrichment > threshold
            & contacts$promoter_contact == TRUE]
    return(distal)
}

distal_contacts <- distal_contact_creation(gr_contacts, promoters)

save(distal_contacts, file = "contact_Rda/distal_contacts.Rda")

#-------------------------Expression Data Processing---------------------------#
#------------------------------------------------------------------------------#


x <- read.table("~/Human_Datasets/H9/Expression/GSM2061386_H9.genes.fpkm.txt",
  header = T)
genes <- import("~/Human_Datasets/Homo_sapiens.GRCh37.75.gtf.gz")
chain <- import.chain("/home/jw18713/Human_Datasets/hg19ToHg38.over.chain")

pos <- match(x$gene_id, genes$gene_id)
x$strand <- NA
x$strand[which(!is.na(pos))] <- strand(genes[pos[!is.na(pos)]])

x$locus <- as.character(x$locus)
chrs <- unlist(strsplit(x$locus, ":"))[c(T,F)]
ranges <- unlist(strsplit(x$locus, ":"))[c(F,T)]
starts <- as.numeric(unlist(strsplit(ranges, "-"))[c(T,F)])
ends <- as.numeric(unlist(strsplit(ranges, "-"))[c(F,T)])

expression_data <- GRanges(seqnames = chrs, ranges=IRanges(starts, ends),
strand = x$strand)
seqlevelsStyle(expression_data) <- "UCSC"
expression_data$Gene_short_name <- as.character(x$gene_short_name)
expression_data$FPKM <- x$FPKM

expression_data <- unlist(liftOver(expression_data, chain))

save(expression_data, file = "expression_data.Rda")
