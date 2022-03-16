library(GenomicRanges)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(PWMEnrich)
library(PWMEnrich.Dmelanogaster.background)
library(MotifDb)

BG3_pred_enh <- get(load("/home/jw18713/Archive_Year1/grown_predicted_enhancers_BG3.Rda")) # BG3_pred_enh
BG3_starr <- get(load("/home/jw18713/Archive_Year1/starr.Rda")) # as BG3_starr

BG3_common_enhancers <- BG3_pred_enh[BG3_pred_enh %over% BG3_starr]
BG3_putative_enhancers <- BG3_pred_enh[!BG3_pred_enh %over% BG3_starr]
BG3_starr_only <- BG3_starr[!BG3_starr %over% BG3_pred_enh]
BG3_putative_enhancers <- BG3_putative_enhancers[width(BG3_putative_enhancers) >= 50]

BG3_common_seq <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, BG3_common_enhancers)
BG3_putative_seq <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, BG3_putative_enhancers)

data(PWMLogn.dm3.MotifDb.Dmel)

registerCoresPWMEnrich(24)

# Super Enhancers, Proximal, Distal Only

if(file.exists("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/BG3_putative_enhancer_PWM_enrichment.RData")){
  load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/BG3_putative_enhancer_PWM_enrichment.RData")
} else{
  BG3_putative_enrichment <- motifEnrichment(BG3_putative_seq, PWMLogn.dm3.MotifDb.Dmel)
  save(BG3_putative_enrichment, file="/home/jw18713/Project1/Paper_Rda_Stuff/BG3/BG3_putative_enhancer_PWM_enrichment.RData")
}

if(file.exists("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/BG3_common_PWM_enrichment.RData")){
  load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/BG3_common_PWM_enrichment.RData")
} else{
  BG3_common_enrichment <- motifEnrichment(BG3_common_seq, PWMLogn.dm3.MotifDb.Dmel)
  save(BG3_common_enrichment, file="/home/jw18713/Project1/Paper_Rda_Stuff/BG3/BG3_common_PWM_enrichment.RData")
}

BG3_putative_report <- groupReport(BG3_putative_enrichment)
BG3_common_report <- groupReport(BG3_common_enrichment)

BG3_putative_report_pvalue05 <- BG3_putative_report[BG3_putative_report$p.value < 0.05]
BG3_common_report_pvalue05 <- BG3_common_report[BG3_common_report$p.value < 0.05]

save(BG3_putative_report_pvalue05, file = "BG3_put_full_report_0.5.Rda")
save(BG3_common_report_pvalue05, file = "BG3_common_full_report_0.5.Rda")

pdf("/home/jw18713/Project1/Paper_Plots/PWM/BG3/BG3_putative_enrichment_report.pdf", width=8, height=11)
plot(BG3_putative_report_pvalue05[1:10], fontsize=7, id.fontsize=5)
dev.off()

pdf("/home/jw18713/Project1/Paper_Plots/PWM/BG3/BG3_common_enrichment_report.pdf", width=8, height=11)
plot(BG3_common_report_pvalue05[1:10], fontsize=7, id.fontsize=5)
dev.off()

write.table(unique(BG3_putative_report_pvalue05$target),
  file="/home/jw18713/Project1/Paper_Plots/PWM/BG3/BG3_putative_report_pvalue05.txt",
  row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(unique(BG3_common_report_pvalue05$target),
  file="/home/jw18713/Project1/Paper_Plots/PWM/BG3/BG3_common_report_pvalue05.txt",
  row.names=FALSE, col.names=FALSE, quote=FALSE)

# get the differential enrichment
BG3_putative_over_common_enrichement_diff <- motifDiffEnrichment(BG3_putative_seq, BG3_common_seq, PWMLogn.dm3.MotifDb.Dmel, verbose=FALSE)
# motifs differentially enriched in the first sequence (with lognormal background correction)
BG3_putative_diff_enrichement <- head(sort(BG3_putative_over_common_enrichement_diff$group.bg, decreasing=TRUE), 100)
BG3_putative_diff_enrichement_names <- names(BG3_putative_diff_enrichement)
ids <- match(names(BG3_putative_diff_enrichement), names(PWMLogn.dm3.MotifDb.Dmel$pwms))
for(i in 1:length(ids)){
  BG3_putative_diff_enrichement_names[i] <- PWMLogn.dm3.MotifDb.Dmel$pwms[[ids[i]]]$name
}

BG3_common_over_putative_enrichement_diff <- motifDiffEnrichment(BG3_common_seq, BG3_putative_seq, PWMLogn.dm3.MotifDb.Dmel, verbose=FALSE)
# motifs differentially enriched in the first sequence (with lognormal background correction)
BG3_common_diff_enrichement <- head(sort(BG3_common_over_putative_enrichement_diff$group.bg, decreasing=TRUE), 100)
BG3_common_diff_enrichement_names <- names(BG3_common_diff_enrichement)
ids <- match(names(BG3_common_diff_enrichement), names(PWMLogn.dm3.MotifDb.Dmel$pwms))
for(i in 1:length(ids)){
BG3_common_diff_enrichement_names[i] <- PWMLogn.dm3.MotifDb.Dmel$pwms[[ids[i]]]$name
}

write(unique(BG3_putative_diff_enrichement_names), file="/home/jw18713/Project1/Paper_Tables/BG3_putative_diff_enriched_TFs.dat")
write(unique(BG3_common_diff_enrichement_names), file="/home/jw18713/Project1/Paper_Tables/BG3_common_diff_enriched_TFs.dat")
