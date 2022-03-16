
library(seqLogo)
library(GenomicRanges)
library(BSgenome.Dmelanogaster.UCSC.dm6)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
library(PWMEnrich)
library(PWMEnrich.Dmelanogaster.background)
library(MotifDb)
library(gridExtra)
# Load the required packages
require(ggplot2)
require(ggseqlogo)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


shared_TFs <- function(x,y){
  xnames <- read.table(x, stringsAsFactors = F)[,1]
  ynames <- read.table(y, stringsAsFactors = F)[,1]
  shared <- xnames[xnames %in% ynames]
}

put_names <- shared_TFs("/home/jw18713/Project1/Paper_Tables/BG3_putative_diff_enriched_TFs.dat",
"/home/jw18713/Project1/Paper_Tables/S2_putative_diff_enriched_TFs.dat")

common_names <- shared_TFs("/home/jw18713/Project1/Paper_Tables/BG3_common_diff_enriched_TFs.dat",
"/home/jw18713/Project1/Paper_Tables/S2_common_diff_enriched_TFs.dat")

put_specific <- put_names[!put_names %in% common_names]
common_specific <- common_names[!common_names %in% put_names]

BG3_put_report <- get(load("BG3_put_full_report_0.5.Rda"))
BG3_common_report <- get(load("BG3_common_full_report_0.5.Rda"))
S2_put_report <- get(load("S2_putative_report.Rda"))
S2_common_report <- get(load("S2_common_report.Rda"))



BG3_put_to_plot <- BG3_put_report[BG3_put_report$target %in% put_specific]
BG3_common_to_plot <- BG3_common_report[BG3_common_report$target %in% common_specific]
S2_put_to_plot <- S2_put_report[S2_put_report$target %in% put_specific]
S2_common_to_plot <- S2_common_report[S2_common_report$target %in% common_specific]





# Next time create these lists by comparing enriched motifs in R...
# this was painful and you should have known better. It isn't a fun
# word search

# This is how you find the right thing from the report
# SE_report_pvalue05[SE_report_pvalue05$target == "ftz-f1"]
# Then take the id to make this list

sciNotationPowerReduced <- function(x, digits = 1, prefix = "") {
  if (length(x) > 1) {
    return(append(sciNotationPower(x[1]), sciNotationPower(x[-1])))
  }
  if (!x) return(0)
  exponent <- floor(log10(x))
  base <- round(x / 10^exponent, digits)
  as.expression(substitute(prefix ~ 10^exponent,
                           list(exponent = exponent, prefix = prefix)))
}

pdf(file = "/home/jw18713/Project1/Paper_Plots/PWM/FigureS6A.pdf",
width = 8, height = 12, pointsize = 14)
par(mar=c(6, 8, 3, 7)+0.1)
par(xpd = NA)

BG3_values <- BG3_put_to_plot$p.value
names(BG3_values) <- BG3_put_to_plot$target
#values <- values[match(unique(as.character(BG3_specific_strong_mRNA[,2])), names(values))]
BG3_values <- BG3_values[!duplicated(names(BG3_values))]
BG3_values <- BG3_values[order(BG3_values, decreasing=F)]
BG3_values <- rev(-log10(BG3_values))

S2_values <- S2_put_to_plot$p.value
names(S2_values) <- S2_put_to_plot$target
#values <- values[match(unique(as.character(S2_specific_strong_mRNA[,2])), names(values))]
S2_values <- S2_values[!duplicated(names(S2_values))]
S2_values <- S2_values[order(S2_values, decreasing=F)]
S2_values <- rev(-log10(S2_values))

combined <- rbind(as.data.frame(t(S2_values)), as.data.frame(t(BG3_values)))
ord_comb <- combined[,order(sapply(data.frame(combined), FUN = mean), decreasing = F)]
colnames(ord_comb)[27] <- "E(spl)my-HLH"
colnames(ord_comb)[23] <- "FoxL1"
colnames(ord_comb)[20] <- "FoxP"
colnames(ord_comb)[17] <- "br"


m <- barplot(as.matrix(ord_comb), col=c(cbbPalette[7:6]), xlim = c(0,200), names.arg = colnames(ord_comb),
             xaxt="n", main="Enriched TFs at Putative Enhancers", xlab="", horiz=TRUE, beside=T, border=NA, las = 1)

#axis(2, labels=colnames(ord_comb), las=2, cex.axis=0.9, tick = FALSE)
axis(BELOW<-1, at=(c(1,10,25,50,100,150,200)), labels=sapply(10^(-c(1,10,25,50,100,150,200)),sciNotationPowerReduced), las=2, cex.axis=1.2, tick = TRUE)
title(xlab="p-value", line=4.5, cex.lab=1.2)
legend(210, 60, legend = c("BG3", "S2"), fill = cbbPalette[7:6], bty="n")

dev.off()

put_ord_comb <- colnames(ord_comb)

pdf(file = "/home/jw18713/Project1/Paper_Plots/PWM/FigureS6B.pdf",
width = 8, height = 12, pointsize = 14)
par(mar=c(6, 8, 3, 7)+0.1)
par(xpd = NA)

BG3_values <- BG3_common_to_plot$p.value
names(BG3_values) <- BG3_common_to_plot$target
#values <- values[match(unique(as.character(BG3_specific_strong_mRNA[,2])), names(values))]
BG3_values <- BG3_values[!duplicated(names(BG3_values))]
BG3_values <- BG3_values[order(BG3_values, decreasing=F)]
BG3_values <- rev(-log10(BG3_values))

S2_values <- S2_common_to_plot$p.value
names(S2_values) <- S2_common_to_plot$target
#values <- values[match(unique(as.character(S2_specific_strong_mRNA[,2])), names(values))]
S2_values <- S2_values[!duplicated(names(S2_values))]
S2_values <- S2_values[order(S2_values, decreasing=F)]
S2_values <- rev(-log10(S2_values))

combined <- rbind(as.data.frame(t(S2_values)), as.data.frame(t(BG3_values)))
ord_comb <- combined[,order(sapply(data.frame(combined), FUN = mean), decreasing = F)]
colnames(ord_comb)[29] <- "dati"
colnames(ord_comb)[26] <- "lov"
colnames(ord_comb)[7] <- "SREBP"
colnames(ord_comb)[5] <- "lrbp18"
colnames(ord_comb)[1] <- "mio"

m <- barplot(as.matrix(ord_comb), col=c(cbbPalette[7:6]), xlim = c(0,150), names.arg = colnames(ord_comb),
             xaxt="n", main="Enriched TFs at Common Enhancers", xlab="", horiz=TRUE, beside=T, border=NA, las = 1)

#axis(2, labels=colnames(ord_comb), las=2, cex.axis=0.9, tick = FALSE)
axis(BELOW<-1, at=(c(1,10,25,50,100,150)), labels=sapply(10^(-c(1,10,25,50,100,150)),sciNotationPowerReduced), las=2, cex.axis=1.2, tick = TRUE)
title(xlab="p-value", line=4.5, cex.lab=1.2)
legend(160, 60, legend = c("BG3", "S2"), fill = cbbPalette[7:6], bty="n")

dev.off()

common_ord_comb <- colnames(ord_comb)

write.table(put_ord_comb, file = "putative_TFs.txt", sep = "\t", row.names = F, quote = F, col.names = F)
write.table(common_ord_comb, file = "common_TFs.txt", sep = "\t", row.names = F, quote = F, col.names = F)







# Everything down here is old stuff
motif_ensmallening <- function(x, report, pval){
  buffer <- rep(0, length(x))
  for (i in seq_along(x)){
    buffer[i] <- report$id[report$target == x[i] & report$p.value <= pval][1]
  }
  return(buffer)
}


put_shared <- motif_ensmallening(put_specific, put_report, 0.05)
common_shared <- motif_ensmallening(common_specific, common_report, 0.05)


pwm_grabber <- function(x){
  pwms <- list()
  for (i in x){
    pwm <- as.list(query(MotifDb,i)[1])
    pwms <- c(pwm, pwms)
  }
  return(pwms)
}

put_pwms <- pwm_grabber(put_shared)
names(put_pwms) <- put_names
common_pwms <- pwm_grabber(common_shared)
names(common_pwms) <- common_names

pdf("/home/jw18713/Project1/Paper_Plots/PWM/std_pwms.pdf",
  width = 14, height = 2, pointsize = 14)
par(cex = 1.2)
ggseqlogo(std_pwms, ncol=3)
dev.off()

pdf("/home/jw18713/Project1/Paper_Plots/PWM/SE_pwms.pdf",
width = 14, height = 6, pointsize = 14)
par(cex = 1.2)
ggseqlogo(SE_pwms, ncol=3)
dev.off()
