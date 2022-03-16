
library(seqLogo)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(PWMEnrich)
library(PWMEnrich.Hsapiens.background)
library(MotifDb)
library(gridExtra)
require(ggplot2)
require(ggseqlogo)
library(genomation)
library(rtracklayer)
library(VennDiagram)


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

readTFs <- function(x){
  names <- read.table(x, stringsAsFactors = F)[,1]
  return(names)
}

putative <-
  get(load("~/Human_Datasets/H9/PWMEnrich100/H9_common_full_report_0.001.Rda"))
common <-
  get(load("~/Human_Datasets/H9/PWMEnrich100/H9_put_full_report_0.001.Rda"))

putative <- putative[!putative$target %in% common$target]
common <- common[!common$target %in% putative$target]

p_top10 <- sort(putative$raw.score, decreasing = T)[1:11]
putative <- putative[putative$raw.score %in% p_top10]

c_top10 <- sort(common$raw.score, decreasing = T)[1:11]
common <- common[common$raw.score %in% c_top10]

# some of the names are not the official symbols
# put_motifs[2] is dl_2, correct symbol for plot is dl
# put_motifs[3] is CG14962, correct symbol for plot is Asciz
# common_motifs[1] CNC::maf-S is CNC with maf-s, id wise I should probably leave it this way
# I'll use the original target names to get the correct motifs though

putative_plot_names <- putative$target
common_plot_names <- common$target

report_id_extract <- function(report, names){
  ids <- rep(0, length(names))
  for (i in seq_along(names)){
    ids[i] <- report$id[report$target == names[i]][1]
  }
  return(ids)
}

put_ids <- report_id_extract(putative, putative_plot_names)
common_ids <- report_id_extract(common, common_plot_names)

pwm_grabber <- function(x){
  pwms <- list()
  for (i in x){
    pwm <- as.list(query(MotifDb,i)[1])
    pwms <- c(pwm, pwms)
  }
  return(pwms)
}

put_pwms <- unique(pwm_grabber(put_ids))
names(put_pwms) <- unique(putative_plot_names)
common_pwms <- unique(pwm_grabber(common_ids))
names(common_pwms) <- unique(common_plot_names)

# Removing Unknown Motifs
put_pwms <- put_pwms[-c(4,10)]
common_pwms <- common_pwms[-6]

pdf("~/Human_Datasets/Plots/Figure5/put_pwms_no_UW.pdf",
  width = 8, height = 8, pointsize = 14)
par(cex = 1.8)
ggseqlogo(put_pwms, ncol=2)
dev.off()

pdf("~/Human_Datasets/Plots/Figure5/common_pwms_no_UW.pdf",
width = 8, height = 10, pointsize = 14)
par(cex = 1.8)
ggseqlogo(common_pwms, ncol=2)
dev.off()
