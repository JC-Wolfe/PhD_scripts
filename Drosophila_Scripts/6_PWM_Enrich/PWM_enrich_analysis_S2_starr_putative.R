library(GenomicRanges)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(PWMEnrich)
library(PWMEnrich.Dmelanogaster.background)
library(MotifDb)

S2_pred_enh <- get(load("/home/jw18713/Archive_Year1/grown_predicted_enhancers_S2.Rda")) # S2_pred_enh
S2_starr <- get(load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/starr.Rda")) # as S2_starr

S2_shared_enhancers <- S2_pred_enh[S2_pred_enh %over% S2_starr]
S2_putative_enhancers <- S2_pred_enh[!S2_pred_enh %over% S2_starr]
S2_starr_only <- S2_starr[!S2_starr %over% S2_pred_enh]
S2_putative_enhancers <- S2_putative_enhancers[width(S2_putative_enhancers) >= 50]

S2_starr_seq <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, S2_starr_only)
S2_putative_seq <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, S2_putative_enhancers)

data(PWMLogn.dm3.MotifDb.Dmel)

registerCoresPWMEnrich(24)

# Super Enhancers, Proximal, Distal Only

if(file.exists("/home/jw18713/Project1/Paper_Rda_Stuff/S2/S2_putative_enhancer_PWM_enrichment.RData")){
  load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/S2_putative_enhancer_PWM_enrichment.RData")
} else{
  S2_putative_enrichment <- motifEnrichment(S2_putative_seq, PWMLogn.dm3.MotifDb.Dmel)
  save(S2_putative_enrichment, file="/home/jw18713/Project1/Paper_Rda_Stuff/S2/S2_putative_enhancer_PWM_enrichment.RData")
}

if(file.exists("/home/jw18713/Project1/Paper_Rda_Stuff/S2/S2_starr_only_PWM_enrichment.RData")){
  load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/S2_starr_only_PWM_enrichment.RData")
} else{
  S2_starr_only_enrichment <- motifEnrichment(S2_starr_seq, PWMLogn.dm3.MotifDb.Dmel)
  save(S2_starr_only_enrichment, file="/home/jw18713/Project1/Paper_Rda_Stuff/S2/S2_starr_only_PWM_enrichment.RData")
}

S2_putative_report <- groupReport(S2_putative_enrichment)
S2_starr_report <- groupReport(S2_starr_only_enrichment)

S2_putative_report_pvalue05 <- S2_putative_report[S2_putative_report$p.value < 0.05]
S2_starr_report_pvalue05 <- S2_starr_report[S2_starr_report$p.value < 0.05]

pdf("/home/jw18713/Project1/Paper_Plots/PWM/S2/S2_putative_enrichment_report.pdf", width=8, height=11)
plot(S2_putative_report_pvalue05[1:10], fontsize=7, id.fontsize=5)
dev.off()

pdf("/home/jw18713/Project1/Paper_Plots/PWM/S2/S2_starr_enrichment_report.pdf", width=8, height=11)
plot(S2_starr_report_pvalue05[1:10], fontsize=7, id.fontsize=5)
dev.off()

write.table(unique(S2_putative_report_pvalue05$target),
  file="/home/jw18713/Project1/Paper_Plots/PWM/S2/SE_proximal_report_pvalue05.txt",
  row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(unique(S2_starr_report_pvalue05$target),
  file="/home/jw18713/Project1/Paper_Plots/PWM/S2/SE_distal_report_pvalue05.txt",
  row.names=FALSE, col.names=FALSE, quote=FALSE)


# get the differential enrichment
S2_putative_over_starr_enrichement_diff <- motifDiffEnrichment(S2_putative_seq, S2_starr_seq, PWMLogn.dm3.MotifDb.Dmel, verbose=FALSE)
# motifs differentially enriched in the first sequence (with lognormal background correction)
S2_putative_diff_enrichement <- head(sort(S2_putative_over_starr_enrichement_diff$group.bg, decreasing=TRUE), 100)
S2_putative_diff_enrichement_names <- names(S2_putative_diff_enrichement)
ids <- match(names(S2_putative_diff_enrichement), names(PWMLogn.dm3.MotifDb.Dmel$pwms))
for(i in 1:length(ids)){
  S2_putative_diff_enrichement_names[i] <- PWMLogn.dm3.MotifDb.Dmel$pwms[[ids[i]]]$name
}

S2_starr_over_putative_enrichement_diff <- motifDiffEnrichment(S2_starr_seq, S2_putative_seq, PWMLogn.dm3.MotifDb.Dmel, verbose=FALSE)
# motifs differentially enriched in the first sequence (with lognormal background correction)
S2_starr_diff_enrichement <- head(sort(S2_starr_over_putative_enrichement_diff$group.bg, decreasing=TRUE), 100)
S2_starr_diff_enrichement_names <- names(S2_starr_diff_enrichement)
ids <- match(names(S2_starr_diff_enrichement), names(PWMLogn.dm3.MotifDb.Dmel$pwms))
for(i in 1:length(ids)){
S2_starr_diff_enrichement_names[i] <- PWMLogn.dm3.MotifDb.Dmel$pwms[[ids[i]]]$name
}


write(unique(S2_putative_diff_enrichement_names), file="/home/jw18713/Project1/Paper_Tables/S2_proximal_diff_enriched_TFs.dat")
write(unique(S2_starr_diff_enrichement_names), file="/home/jw18713/Project1/Paper_Tables/S2_distal_only_diff_enriched_TFs.dat")
