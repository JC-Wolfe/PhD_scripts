library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(PWMEnrich)
# Importing libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("PWMEnrich.Hsapiens.background")
library(PWMEnrich.Hsapiens.background)
library(MotifDb)

load("~/Human_Datasets/H9/Processed_datasets/H9_full(widths).Rda")

H9_putative <- reduce(H9_full[H9_full$classification == "Putative"])
H9_common <- reduce(H9_full[H9_full$classification == "Common"])

H9_putative <- H9_putative[width(H9_putative) >= 200]
H9_common <- H9_common[width(H9_common) >= 200]

H9_putative_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, H9_putative)
H9_common_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, H9_common)

data(PWMLogn.hg19.MotifDb.Hsap)

registerCoresPWMEnrich(40)

# Super Enhancers, Proximal, Distal Only
setwd("~/Human_Datasets/H9/PWMEnrich100/")

if(file.exists("H9_putative_enhancer_PWM_enrichment.RData")){
  load("H9_putative_enhancer_PWM_enrichment.RData")
} else{
  H9_putative_enrichment <- motifEnrichment(H9_putative_seq, PWMLogn.hg19.MotifDb.Hsap)
  save(H9_putative_enrichment, file="H9_putative_enhancer_PWM_enrichment.RData")
}

if(file.exists("H9_common_PWM_enrichment.RData")){
  load("H9_common_PWM_enrichment.RData")
} else{
  H9_common_enrichment <- motifEnrichment(H9_common_seq, PWMLogn.hg19.MotifDb.Hsap)
  save(H9_common_enrichment, file="H9_common_PWM_enrichment.RData")
}

H9_putative_report <- groupReport(H9_putative_enrichment)
H9_common_report <- groupReport(H9_common_enrichment)

H9_putative_report_pvalue001 <- H9_putative_report[H9_putative_report$p.value <= 0.001]
H9_common_report_pvalue001 <- H9_common_report[H9_common_report$p.value <= 0.001]

save(H9_putative_report_pvalue001, file = "H9_put_full_report_0.001.Rda")
save(H9_common_report_pvalue001, file = "H9_common_full_report_0.001.Rda")

pdf("H9_putative_enrichment_report.pdf", width=8, height=11)
plot(H9_putative_report_pvalue001[1:10], fontsize=7, id.fontsize=5)
dev.off()

pdf("H9_common_enrichment_report.pdf", width=8, height=11)
plot(H9_common_report_pvalue001[1:10], fontsize=7, id.fontsize=5)
dev.off()

write.table(unique(H9_putative_report_pvalue001$target),
  file="H9_putative_report_pvalue001.txt",
  row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(unique(H9_common_report_pvalue001$target),
  file="H9_common_report_pvalue001.txt",
  row.names=FALSE, col.names=FALSE, quote=FALSE)

# get the differential enrichment
H9_putative_over_common_enrichement_diff <- motifDiffEnrichment(H9_putative_seq, H9_common_seq, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE)
# motifs differentially enriched in the first sequence (with lognormal background correction)
H9_putative_diff_enrichement <- head(sort(H9_putative_over_common_enrichement_diff$group.bg, decreasing=TRUE), 100)
H9_putative_diff_enrichement_names <- names(H9_putative_diff_enrichement)
ids <- match(names(H9_putative_diff_enrichement), names(PWMLogn.hg19.MotifDb.Hsap$pwms))
for(i in 1:length(ids)){
  H9_putative_diff_enrichement_names[i] <- PWMLogn.hg19.MotifDb.Hsap$pwms[[ids[i]]]$name
}

H9_common_over_putative_enrichement_diff <- motifDiffEnrichment(H9_common_seq, H9_putative_seq, PWMLogn.hg19.MotifDb.Hsap, verbose=FALSE)
# motifs differentially enriched in the first sequence (with lognormal background correction)
H9_common_diff_enrichement <- head(sort(H9_common_over_putative_enrichement_diff$group.bg, decreasing=TRUE), 100)
H9_common_diff_enrichement_names <- names(H9_common_diff_enrichement)
ids <- match(names(H9_common_diff_enrichement), names(PWMLogn.hg19.MotifDb.Hsap$pwms))
for(i in 1:length(ids)){
H9_common_diff_enrichement_names[i] <- PWMLogn.hg19.MotifDb.Hsap$pwms[[ids[i]]]$name
}

write(unique(H9_putative_diff_enrichement_names), file="H9_putative_diff_enriched_TFs.dat")
write(unique(H9_common_diff_enrichement_names), file="H9_common_diff_enriched_TFs.dat")
