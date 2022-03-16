
library(seqLogo)
library(GenomicRanges)
library(BSgenome.Dmelanogaster.UCSC.dm6)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("PWMEnrich")
library(PWMEnrich)
library(PWMEnrich.Dmelanogaster.background)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MotifDb")
library(MotifDb)
library(gridExtra)
# Load the required packages
require(ggplot2)
require(ggseqlogo)

shared_TFs <- function(x,y){
  xnames <- read.table(x, stringsAsFactors = F)[,1]
  ynames <- read.table(y, stringsAsFactors = F)[,1]
  shared <- xnames[xnames %in% ynames]
}

std_names <- shared_TFs("/home/jw18713/Project1/Paper_Plots/PWM/BG3_std_specific.txt",
"/home/jw18713/Project1/Paper_Plots/PWM/S2_std_specific.txt")

SE_names <- shared_TFs("/home/jw18713/Project1/Paper_Plots/PWM/BG3_SE_specific.txt",
"/home/jw18713/Project1/Paper_Plots/PWM/S2_SE_specific.txt")

# Next time create these lists by comparing enriched motifs in R...
# this was painful and you should have known better. It isn't a fun
# word search

# This is how you find the right thing from the report
# SE_report_pvalue05[SE_report_pvalue05$target == "ftz-f1"]
# Then take the id to make this list

std_shared <- c(
"ftz-f1_FlyReg_FBgn0001078",
"M5177_1.02",
"Atf6_SANGER_5_FBgn0033010"
)

SE_shared <- c(
"BtbVII_SANGER_5_FBgn0263108",
"jigr1_SANGER_5_FBgn0039350",
"Cf2_II",
"Blimp-1_NAR_FBgn0035625",
"Sox14_SANGER_10_FBgn0005612",
"eg_SANGER_5_FBgn0000560",
"CNC::maf-S",
"Atf-2_SANGER_5_FBgn0050420",
"CG6276_SANGER_5_FBgn0038316"
)

pwm_grabber <- function(x){
  pwms <- list()
  for (i in x){
    pwm <- as.list(query(MotifDb,i)[1])
    pwms <- c(pwm, pwms)
  }
  return(pwms)
}

std_pwms <- pwm_grabber(std_shared)
names(std_pwms) <- std_names
SE_pwms <- pwm_grabber(SE_shared)
names(SE_pwms) <- SE_names

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
