
library(GenomicRanges)
library(rtracklayer)

load("/home/jw18713/rerun_glists/starr/BG3_starr_dist.Rda")
load("/home/jw18713/rerun_glists/starr/BG3_starr_prox.Rda")
load("/home/jw18713/rerun_glists/bg/BG3_bg_prox.Rda")
load("/home/jw18713/rerun_glists/bg/BG3_bg_dist.Rda")


load("/home/jw18713/rerun_glists/starr/S2_starr_dist.Rda")
load("/home/jw18713/rerun_glists/starr/S2_starr_prox.Rda")
load("/home/jw18713/rerun_glists/bg/S2_bg_prox.Rda")
load("/home/jw18713/rerun_glists/bg/S2_bg_dist.Rda")

ann_genome <- get(load(file = "/home/jw18713/Archive_Year1/annotated_GRange_FI_done.Rda"))
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Proximal BG3 (BG3P)

# Now which of these don't have any expresion?
no_expression_BG3P <- BG3_starr_prox[lengths(BG3_starr_prox) == 0]

# And which of these are expressed?
enhancers_contacted_BG3P <- BG3_starr_prox[lengths(BG3_starr_prox) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$FPKM, na.rm = T))
}

max_to_plot_BG3P <- sapply(enhancers_contacted_BG3P, maxExp)

no_expression_max_BG3P <- c(no_expression_BG3P, enhancers_contacted_BG3P[max_to_plot_BG3P == 0])
expression_max_BG3P <- max_to_plot_BG3P[max_to_plot_BG3P > 0]

# DISTAL ONLY BG3 (BG3D)

# Now which of these don't have any expresion?
no_expression_BG3D <- BG3_starr_dist[lengths(BG3_starr_dist) == 0]

# And which of these are expressed?
enhancers_contacted_BG3D <- BG3_starr_dist[lengths(BG3_starr_dist) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$FPKM, na.rm = T))
}

max_to_plot_BG3D <- sapply(enhancers_contacted_BG3D, maxExp)

no_expression_max_BG3D <- c(no_expression_BG3D, enhancers_contacted_BG3D[max_to_plot_BG3D == 0])
expression_max_BG3D <- max_to_plot_BG3D[max_to_plot_BG3D > 0]

# BG3 Proximal background (BG3_bg_p)

# Now which of these don't have any expresion?
no_expression_BG3_bg_p <- BG3_bg_prox[lengths(BG3_bg_prox) == 0]

# And which of these are expressed?
enhancers_contacted_BG3_bg_p <- BG3_bg_prox[lengths(BG3_bg_prox) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$FPKM, na.rm = T))
}

max_to_plot_BG3_bg_p <- sapply(enhancers_contacted_BG3_bg_p, maxExp)

no_expression_max_BG3_bg_p <- c(no_expression_BG3_bg_p, enhancers_contacted_BG3_bg_p[max_to_plot_BG3_bg_p == 0])
expression_max_BG3_bg_p <- max_to_plot_BG3_bg_p[max_to_plot_BG3_bg_p > 0]



# BG3 DISTAL ONLY backgound (BG3_bg_do)

# Now which of these don't have any expresion?
no_expression_BG3_bg_do <- BG3_bg_dist[lengths(BG3_bg_dist) == 0]

# And which of these are expressed?
enhancers_contacted_BG3_bg_do <- BG3_bg_dist[lengths(BG3_bg_dist) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$FPKM, na.rm = T))
}

max_to_plot_BG3_bg_do <- sapply(enhancers_contacted_BG3_bg_do, maxExp)

no_expression_max_BG3_bg_do <- c(no_expression_BG3_bg_do, enhancers_contacted_BG3_bg_do[max_to_plot_BG3_bg_do == 0])
expression_max_BG3_bg_do <- max_to_plot_BG3_bg_do[max_to_plot_BG3_bg_do > 0]




BG3_em <- matrix(0, 4, 2, dimnames = list(c("Distal Only", "Background\nDistal Only", "Proximal", "Background\nProximal"),
  c("Expressed", "Not Expressed")))

BG3_em[1,1] <- length(expression_max_BG3D)
BG3_em[1,2] <- length(no_expression_max_BG3D)
BG3_em[2,1] <- length(expression_max_BG3_bg_do)
BG3_em[2,2] <- length(no_expression_max_BG3_bg_do)
BG3_em[3,1] <- length(expression_max_BG3P)
BG3_em[3,2] <- length(no_expression_max_BG3P)
BG3_em[4,1] <- length(expression_max_BG3_bg_p)
BG3_em[4,2] <- length(no_expression_max_BG3_bg_p)

starr_dist <- BG3_em[1,]
bg_dist <- BG3_em[2,]
starr_prox <- BG3_em[3,]
bg_prox <- BG3_em[4,]


fisher.test(rbind(starr_prox, starr_dist))
fisher.test(rbind(starr_prox, bg_prox))
fisher.test(rbind(starr_prox, bg_dist))
fisher.test(rbind(starr_dist, bg_prox))
fisher.test(rbind(starr_dist, bg_dist))
fisher.test(rbind(bg_prox, bg_dist))



# Proximal S2 (S2P)

# Now which of these don't have any expresion?
no_expression_S2P <- S2_starr_prox[lengths(S2_starr_prox) == 0]

# And which of these are expressed?
enhancers_contacted_S2P <- S2_starr_prox[lengths(S2_starr_prox) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$FPKM, na.rm = T))
}

max_to_plot_S2P <- sapply(enhancers_contacted_S2P, maxExp)

no_expression_max_S2P <- c(no_expression_S2P, enhancers_contacted_S2P[max_to_plot_S2P == 0])
expression_max_S2P <- max_to_plot_S2P[max_to_plot_S2P > 0]

# DISTAL ONLY S2 (S2D)

# Now which of these don't have any expresion?
no_expression_S2D <- S2_starr_dist[lengths(S2_starr_dist) == 0]

# And which of these are expressed?
enhancers_contacted_S2D <- S2_starr_dist[lengths(S2_starr_dist) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$FPKM, na.rm = T))
}

max_to_plot_S2D <- sapply(enhancers_contacted_S2D, maxExp)

no_expression_max_S2D <- c(no_expression_S2D, enhancers_contacted_S2D[max_to_plot_S2D == 0])
expression_max_S2D <- max_to_plot_S2D[max_to_plot_S2D > 0]

# S2 Proximal background (S2_bg_p)

# Now which of these don't have any expresion?
no_expression_S2_bg_p <- S2_bg_prox[lengths(S2_bg_prox) == 0]

# And which of these are expressed?
enhancers_contacted_S2_bg_p <- S2_bg_prox[lengths(S2_bg_prox) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$FPKM, na.rm = T))
}

max_to_plot_S2_bg_p <- sapply(enhancers_contacted_S2_bg_p, maxExp)

no_expression_max_S2_bg_p <- c(no_expression_S2_bg_p, enhancers_contacted_S2_bg_p[max_to_plot_S2_bg_p == 0])
expression_max_S2_bg_p <- max_to_plot_S2_bg_p[max_to_plot_S2_bg_p > 0]



# S2 DISTAL ONLY backgound (S2_bg_do)

# Now which of these don't have any expresion?
no_expression_S2_bg_do <- S2_bg_dist[lengths(S2_bg_dist) == 0]

# And which of these are expressed?
enhancers_contacted_S2_bg_do <- S2_bg_dist[lengths(S2_bg_dist) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$FPKM, na.rm = T))
}

max_to_plot_S2_bg_do <- sapply(enhancers_contacted_S2_bg_do, maxExp)

no_expression_max_S2_bg_do <- c(no_expression_S2_bg_do, enhancers_contacted_S2_bg_do[max_to_plot_S2_bg_do == 0])
expression_max_S2_bg_do <- max_to_plot_S2_bg_do[max_to_plot_S2_bg_do > 0]




S2_em <- matrix(0, 4, 2, dimnames = list(c("Distal Only", "Background\nDistal Only", "Proximal", "Background\nProximal"),
  c("Expressed", "Not Expressed")))

S2_em[1,1] <- length(expression_max_S2D)
S2_em[1,2] <- length(no_expression_max_S2D)
S2_em[2,1] <- length(expression_max_S2_bg_do)
S2_em[2,2] <- length(no_expression_max_S2_bg_do)
S2_em[3,1] <- length(expression_max_S2P)
S2_em[3,2] <- length(no_expression_max_S2P)
S2_em[4,1] <- length(expression_max_S2_bg_p)
S2_em[4,2] <- length(no_expression_max_S2_bg_p)

starr_dist <- S2_em[1,]
bg_dist <- S2_em[2,]
starr_prox <- S2_em[3,]
bg_prox <- S2_em[4,]



fisher.test(rbind(starr_prox, starr_dist))
fisher.test(rbind(starr_prox, bg_prox))
fisher.test(rbind(starr_dist, bg_dist))


#fisher.test(S2_em)

BG3_per <- BG3_em[,1]
BG3_per[1] <- BG3_em[1,1]/(BG3_em[1,1]+BG3_em[1,2])*100
BG3_per[2] <- BG3_em[2,1]/(BG3_em[2,1]+BG3_em[2,2])*100
BG3_per[3] <- BG3_em[3,1]/(BG3_em[3,1]+BG3_em[3,2])*100
BG3_per[4] <- BG3_em[4,1]/(BG3_em[4,1]+BG3_em[4,2])*100

S2_per <- S2_em[,1]
S2_per[1] <- S2_em[1,1]/(S2_em[1,1]+S2_em[1,2])*100
S2_per[2] <- S2_em[2,1]/(S2_em[2,1]+S2_em[2,2])*100
S2_per[3] <- S2_em[3,1]/(S2_em[3,1]+S2_em[3,2])*100
S2_per[4] <- S2_em[4,1]/(S2_em[4,1]+S2_em[4,2])*100

BG3_plot <- matrix(0,2,2)
colnames(BG3_plot) <- c("zDistal Only", "aProximal")
rownames(BG3_plot) <- c("BG3", "Background")
BG3_plot[1,1] <- BG3_per[1]
BG3_plot[1,2] <- BG3_per[3]
BG3_plot[2,1] <- BG3_per[2]
BG3_plot[2,2] <- BG3_per[4]

S2_plot <- matrix(0,2,2)
colnames(S2_plot) <- c("zDistal Only", "aProximal")
rownames(S2_plot) <- c("S2", "Background")
S2_plot[1,1] <- S2_per[1]
S2_plot[1,2] <- S2_per[3]
S2_plot[2,1] <- S2_per[2]
S2_plot[2,2] <- S2_per[4]

pdf(file = "/home/jw18713/Project1/Paper_Plots/FigureS4/FigureS4_Horiz_BG.pdf",
width = 32, height = 12, pointsize = 14)
par(mar=c(4,5.5,4,15), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)

  barplot(BG3_plot[2:1,2:1], beside=T, horiz = T, col = cbbPalette[c(2,7,3,6)],
      xlim = c(0,100), xlab = "Expressed Contacted Genes", xaxt = "n",
      names.arg = c("",""), main = "BG3 Contacted Regions\nWith Expression")
    axis(1, seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"), las=2)
    legend(110,3.5,c("BG3 Distal Only",
                    "Background Distal Only",
                    "BG3 Proximal",
                    "Background Proximal"), fill=cbbPalette[c(6,3,7,2)], bty='n')
    # starr proximal and proximal bg
      lines(c(max(BG3_plot[,2]) + 1.5, max(BG3_plot[,2]) + 3), c(2.5,2.5))
      lines(c(max(BG3_plot[,2]) + 3, max(BG3_plot[,2]) + 3), c(1.5,2.5))
      lines(c(max(BG3_plot[,2]) + 1.5, max(BG3_plot[,2]) + 3), c(1.5,1.5))
      text(x = max(BG3_plot[,2]) + 4, y = 1.5 + (2.5-1.5)/2, "***", srt = 90)
    # starr distal only and starr proximal
      lines(c(max(BG3_plot[,1]) + 4.5, max(BG3_plot[,1]) + 6), c(5.5,5.5))
      lines(c(max(BG3_plot[,1]) + 6, max(BG3_plot[,1]) + 6), c(2.5,5.5))
      lines(c(max(BG3_plot[,1]) + 4.5, max(BG3_plot[,1]) + 6), c(2.5,2.5))
      text(x = max(BG3_plot[,1]) + 7, y = 2.5 + (5.5-2.5)/2, "***", srt = 90)


  barplot(BG3_plot[2:1,2:1], beside=T, horiz = T, col = cbbPalette[c(2,7,3,6)],
      xlim = c(0,100), xlab = "Expressed Contacted Genes", xaxt = "n",
      names.arg = c("",""), main = "BG3 Contacted Regions\nWith Expression")
    axis(1, seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"), las=2)
    legend(110,3.5,c("BG3 Distal Only",
                    "Background Distal Only",
                    "BG3 Proximal",
                    "Background Proximal"), fill=cbbPalette[c(6,3,7,2)], bty='n')
    # starr proximal and proximal bg
      lines(c(max(S2_plot[,2]) + 1.5, max(S2_plot[,2]) + 3), c(2.5,2.5))
      lines(c(max(S2_plot[,2]) + 3, max(S2_plot[,2]) + 3), c(1.5,2.5))
      lines(c(max(S2_plot[,2]) + 1.5, max(S2_plot[,2]) + 3), c(1.5,1.5))
      text(x = max(S2_plot[,2]) + 4, y = 1.5 + (2.5-1.5)/2, "***", srt = 90)
    # starr distal only and starr proximal
      lines(c(max(S2_plot[,1]) + 4.5, max(S2_plot[,1]) + 6), c(5.5,5.5))
      lines(c(max(S2_plot[,1]) + 6, max(S2_plot[,1]) + 6), c(2.5,5.5))
      lines(c(max(S2_plot[,1]) + 4.5, max(S2_plot[,1]) + 6), c(2.5,2.5))
      text(x = max(S2_plot[,1]) + 7, y = 2.5 + (5.5-2.5)/2, "***", srt = 90)

dev.off()
