library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)

setwd("~/Review_Paper_Plots/rda_stuff")

# load putative enhancers
load("putative.Rda")

# load common enhancers
load("common.Rda")

all_enhancers <- c(putative, common)

# load BG3 Hi-C data as gr_contacts
load("/home/jw18713/Project1/data/BG3_enriched_contacts_GR.Rda")
gr_contacts <- gr_contacts[!start(gr_contacts) == gr_contacts$start &
  !end(gr_contacts) == gr_contacts$start]

contacted <- GRanges(seqnames = gr_contacts$chr, ranges = IRanges(
  start = gr_contacts$start, end = gr_contacts$end))

seqlevelsStyle(contacted) <- "UCSC"

# Solving the mcol name problem from the Hi-C data
gr_fixed <- GRanges(seqnames = seqnames(gr_contacts), ranges = IRanges(
  start = start(gr_contacts), end = end(gr_contacts)))
gr_fixed$chr <- gr_contacts$chr
gr_fixed$start_point <- gr_contacts$start
gr_fixed$end_point <- gr_contacts$end
gr_fixed$enrichment <- gr_contacts$enrichment


# getting putative enhancers that have overlaps
pover <- findOverlaps(putative, gr_fixed)
plist <- split(contacted[subjectHits(pover)], queryHits(pover))

# same for common enhancers
cover <- findOverlaps(common, gr_fixed)
clist <- split(contacted[subjectHits(cover)], queryHits(cover))

# Creating a function to count unique enhancers contacted
enhcounts <- function(x, enhancers){
  o <- findOverlaps(x, enhancers)
  return(length(unique(subjectHits(o))))
}

# Contact counts for putative and common enhancers
p_enh_contacts <- lapply(plist, enhcounts, enhancers = all_enhancers)
c_enh_contacts <- lapply(clist, enhcounts, enhancers = all_enhancers)

# Unlisting for plotting
BG3_p_plot <- unlist(p_enh_contacts)
BG3_c_plot <- unlist(c_enh_contacts)

# Right, this may be lazy but I'm just going to combine them by redoing this for S2
# Should be fine as long as I keep the BG3 ones above
# S2 stuff
load("S2_putative.Rda")
load("S2_common.Rda")
all_enhancers <- c(putative, common)

# S2 overlaps stuff
load("~/Project1/data/S2_enriched_contacts_GR.Rda")
gr_contacts <- gr_contacts[!start(gr_contacts) == gr_contacts$start &
  !end(gr_contacts) == gr_contacts$start]

contacted <- GRanges(seqnames = gr_contacts$chr, ranges = IRanges(
  start = gr_contacts$start, end = gr_contacts$end))

seqlevelsStyle(contacted) <- "UCSC"

# Solving the mcol name problem from the S2 Hi-C data
gr_fixed <- GRanges(seqnames = seqnames(gr_contacts), ranges = IRanges(
  start = start(gr_contacts), end = end(gr_contacts)))
gr_fixed$chr <- gr_contacts$chr
gr_fixed$start_point <- gr_contacts$start
gr_fixed$end_point <- gr_contacts$end
gr_fixed$enrichment <- gr_contacts$enrichment


# getting putative S2 enhancers that have overlaps
pover <- findOverlaps(putative, gr_fixed)
plist <- split(contacted[subjectHits(pover)], queryHits(pover))

# same for common S2 enhancers
cover <- findOverlaps(common, gr_fixed)
clist <- split(contacted[subjectHits(cover)], queryHits(cover))

# Contact counts for putative and common enhancers in S2 cells
p_enh_contacts <- lapply(plist, enhcounts, enhancers = all_enhancers)
c_enh_contacts <- lapply(clist, enhcounts, enhancers = all_enhancers)

# Unlisting for plotting
S2_p_plot <- unlist(p_enh_contacts)
S2_c_plot <- unlist(c_enh_contacts)








setwd("~/Review_Paper_Plots/enh_contacts")

pdf(file = "Enhancer_contacts.pdf",
width = 18, height = 8, pointsize = 14)
par(mar=c(5,5.5,5,6), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)

      # BG3
      # BG3 putative
      hist(BG3_c_plot, main="BG3 Predicted Enhancer Contacts",
        xlab="Predicted Enhancers Contacted", col="#D55E0080", cex.main=2, cex.lab=1.5,
        ylim = c(0,6000), xlim = c(0,1500), breaks = 20, ylab=NA)

      #BG3 common
      hist(BG3_p_plot, col="#0072B280", add= T, breaks = 20)

      # S2
      # S2 putative
      hist(S2_c_plot, main="S2 Predicted Enhancer Contacts",
        xlab="Predicted Enhancers Contacted", col="#D55E0080", cex.main=2, cex.lab=1.5,
        ylim = c(0,7000), xlim = c(0,3500), breaks = 20, ylab=NA)

      # S2 common
      hist(S2_p_plot, col="#0072B280", add= T, breaks = 20)

      legend(4000,4750,c("Common", "Putative"), fill=c("#D55E0080", "#0072B280"), bty="n")

dev.off()


hist_density <- function(x, breakmax, n = 20){
  h <- hist(x, plot = F, breaks = seq(0, breakmax, by = breakmax / n))
  h$density <- h$counts/sum(h$counts)*100
  return(h)
}

# To plot, you'll need to plot(hist_density(x, breakmax = n), freq = F)

pdf(file = "Enhancer_contacts_density.pdf",
width = 18, height = 8, pointsize = 14)
par(mar=c(5,5.5,5,6), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)

      # BG3
      # BG3 putative
      plot(hist_density(BG3_c_plot, breakmax = 2000),
        main="BG3 Predicted Enhancer Contacts", yaxt = "n",
        xlab="Predicted Enhancers Contacted", col="#D55E0080", cex.main=2, cex.lab=1.5,
        ylim = c(0, 80), xlim = c(0,2000), ylab=NA, freq = F)
        axis(2, seq(0, 80, 10), labels = paste0(seq(0, 80, 10), "%"), las=2)

      #BG3 common
      plot(hist_density(BG3_p_plot, breakmax = 2000), col="#0072B280", add= T, freq = F)

      # S2
      # S2 putative
      plot(hist_density(S2_c_plot, breakmax = 4000),
        main="S2 Predicted Enhancer Contacts", yaxt = "n",
        xlab="Predicted Enhancers Contacted", col="#D55E0080", cex.main=2, cex.lab=1.5,
        ylim = c(0, 80), xlim = c(0,4000), ylab=NA, freq = F)
        axis(2, seq(0, 80, 10), labels = paste0(seq(0, 80, 10), "%"), las=2)

      # S2 common
      plot(hist_density(S2_p_plot, breakmax = 4000), col="#0072B280", add= T, freq = F)

      legend(4000,45,c("Common", "Putative"), fill=c("#D55E0080", "#0072B280"), bty="n")

dev.off()
