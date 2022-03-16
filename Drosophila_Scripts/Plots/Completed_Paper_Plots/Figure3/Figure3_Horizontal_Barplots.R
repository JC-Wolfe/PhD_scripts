
library(GenomicRanges)
library(rtracklayer)

load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/BG3_distal_only_glist.Rda")
load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/BG3_proximal_enhancers_glist.Rda")
load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/BG3_proximal_background_glist.Rda")
load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/BG3_background_distal_only_glist.Rda")


load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/S2_distal_only_glist.Rda")
load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/S2_proximal_enhancer_glist.Rda")
load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/S2_proximal_background_glist.Rda")
load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/S2_background_distal_only_glist.Rda")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Proximal BG3 (BG3P)

# Now which of these don't have any expresion?
no_expression_BG3P <- BG3_proximal_enhancers_glist[lengths(BG3_proximal_enhancers_glist) == 0]

# And which of these are expressed?
enhancers_contacted_BG3P <- BG3_proximal_enhancers_glist[lengths(BG3_proximal_enhancers_glist) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$BG3_FPKM, na.rm = T))
}

max_to_plot_BG3P <- sapply(enhancers_contacted_BG3P, maxExp)

no_expression_max_BG3P <- c(no_expression_BG3P, enhancers_contacted_BG3P[max_to_plot_BG3P == 0])
expression_max_BG3P <- max_to_plot_BG3P[max_to_plot_BG3P > 0]

# DISTAL ONLY BG3 (BG3D)

# Now which of these don't have any expresion?
no_expression_BG3D <- BG3_distal_only_glist[lengths(BG3_distal_only_glist) == 0]

# And which of these are expressed?
enhancers_contacted_BG3D <- BG3_distal_only_glist[lengths(BG3_distal_only_glist) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$BG3_FPKM, na.rm = T))
}

max_to_plot_BG3D <- sapply(enhancers_contacted_BG3D, maxExp)

no_expression_max_BG3D <- c(no_expression_BG3D, enhancers_contacted_BG3D[max_to_plot_BG3D == 0])
expression_max_BG3D <- max_to_plot_BG3D[max_to_plot_BG3D > 0]

# BG3 Proximal background (BG3_bg_p)

# Now which of these don't have any expresion?
no_expression_BG3_bg_p <- BG3_proximal_background_glist[lengths(BG3_proximal_background_glist) == 0]

# And which of these are expressed?
enhancers_contacted_BG3_bg_p <- BG3_proximal_background_glist[lengths(BG3_proximal_background_glist) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$BG3_FPKM, na.rm = T))
}

max_to_plot_BG3_bg_p <- sapply(enhancers_contacted_BG3_bg_p, maxExp)

no_expression_max_BG3_bg_p <- c(no_expression_BG3_bg_p, enhancers_contacted_BG3_bg_p[max_to_plot_BG3_bg_p == 0])
expression_max_BG3_bg_p <- max_to_plot_BG3_bg_p[max_to_plot_BG3_bg_p > 0]



# BG3 DISTAL ONLY backgound (BG3_bg_do)

# Now which of these don't have any expresion?
no_expression_BG3_bg_do <- BG3_background_distal_only_glist[lengths(BG3_background_distal_only_glist) == 0]

# And which of these are expressed?
enhancers_contacted_BG3_bg_do <- BG3_background_distal_only_glist[lengths(BG3_background_distal_only_glist) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$BG3_FPKM, na.rm = T))
}

max_to_plot_BG3_bg_do <- sapply(enhancers_contacted_BG3_bg_do, maxExp)

no_expression_max_BG3_bg_do <- c(no_expression_BG3_bg_do, enhancers_contacted_BG3_bg_do[max_to_plot_BG3_bg_do == 0])
expression_max_BG3_bg_do <- max_to_plot_BG3_bg_do[max_to_plot_BG3_bg_do > 0]




BG3_em <- matrix(0, 4, 2, dimnames = list(c("Proximal", "Distal Only", "Background\nProximal", "Background\nDistal Only"),
  c("Expressed", "Not Expressed")))

BG3_em[1,1] <- length(expression_max_BG3P)
BG3_em[1,2] <- length(no_expression_max_BG3P)
BG3_em[2,1] <- length(expression_max_BG3D)
BG3_em[2,2] <- length(no_expression_max_BG3D)
BG3_em[3,1] <- length(expression_max_BG3_bg_p)
BG3_em[3,2] <- length(no_expression_max_BG3_bg_p)
BG3_em[4,1] <- length(expression_max_BG3_bg_do)
BG3_em[4,2] <- length(no_expression_max_BG3_bg_do)

put_prox <- BG3_em[1,]
put_dist <- BG3_em[2,]
bg_prox <- BG3_em[3,]
bg_dist <- BG3_em[4,]


fisher.test(rbind(put_prox, put_dist))
fisher.test(rbind(put_prox, bg_prox))
fisher.test(rbind(put_prox, bg_dist))
fisher.test(rbind(put_dist, bg_prox))
fisher.test(rbind(put_dist, bg_dist))
fisher.test(rbind(bg_prox, bg_dist))



# Proximal S2 (S2P)

# Now which of these don't have any expresion?
no_expression_S2P <- S2_proximal_enhancers_glist[lengths(S2_proximal_enhancers_glist) == 0]

# And which of these are expressed?
enhancers_contacted_S2P <- S2_proximal_enhancers_glist[lengths(S2_proximal_enhancers_glist) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$S2_FPKM, na.rm = T))
}

max_to_plot_S2P <- sapply(enhancers_contacted_S2P, maxExp)

no_expression_max_S2P <- c(no_expression_S2P, enhancers_contacted_S2P[max_to_plot_S2P == 0])
expression_max_S2P <- max_to_plot_S2P[max_to_plot_S2P > 0]

# DISTAL ONLY S2 (S2D)

# Now which of these don't have any expresion?
no_expression_S2D <- S2_distal_only_glist[lengths(S2_distal_only_glist) == 0]

# And which of these are expressed?
enhancers_contacted_S2D <- S2_distal_only_glist[lengths(S2_distal_only_glist) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$S2_FPKM, na.rm = T))
}

max_to_plot_S2D <- sapply(enhancers_contacted_S2D, maxExp)

no_expression_max_S2D <- c(no_expression_S2D, enhancers_contacted_S2D[max_to_plot_S2D == 0])
expression_max_S2D <- max_to_plot_S2D[max_to_plot_S2D > 0]

# S2 Proximal background (S2_bg_p)

# Now which of these don't have any expresion?
no_expression_S2_bg_p <- S2_proximal_background_glist[lengths(S2_proximal_background_glist) == 0]

# And which of these are expressed?
enhancers_contacted_S2_bg_p <- S2_proximal_background_glist[lengths(S2_proximal_background_glist) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$S2_FPKM, na.rm = T))
}

max_to_plot_S2_bg_p <- sapply(enhancers_contacted_S2_bg_p, maxExp)

no_expression_max_S2_bg_p <- c(no_expression_S2_bg_p, enhancers_contacted_S2_bg_p[max_to_plot_S2_bg_p == 0])
expression_max_S2_bg_p <- max_to_plot_S2_bg_p[max_to_plot_S2_bg_p > 0]



# S2 DISTAL ONLY backgound (S2_bg_do)

# Now which of these don't have any expresion?
no_expression_S2_bg_do <- S2_background_distal_only_glist[lengths(S2_background_distal_only_glist) == 0]

# And which of these are expressed?
enhancers_contacted_S2_bg_do <- S2_background_distal_only_glist[lengths(S2_background_distal_only_glist) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$S2_FPKM, na.rm = T))
}

max_to_plot_S2_bg_do <- sapply(enhancers_contacted_S2_bg_do, maxExp)

no_expression_max_S2_bg_do <- c(no_expression_S2_bg_do, enhancers_contacted_S2_bg_do[max_to_plot_S2_bg_do == 0])
expression_max_S2_bg_do <- max_to_plot_S2_bg_do[max_to_plot_S2_bg_do > 0]




S2_em <- matrix(0, 4, 2, dimnames = list(c("Proximal", "Distal Only", "Background\nProximal", "Background\nDistal Only"),
  c("Expressed", "Not Expressed")))

S2_em[1,1] <- length(expression_max_S2P)
S2_em[1,2] <- length(no_expression_max_S2P)
S2_em[2,1] <- length(expression_max_S2D)
S2_em[2,2] <- length(no_expression_max_S2D)
S2_em[3,1] <- length(expression_max_S2_bg_p)
S2_em[3,2] <- length(no_expression_max_S2_bg_p)
S2_em[4,1] <- length(expression_max_S2_bg_do)
S2_em[4,2] <- length(no_expression_max_S2_bg_do)

S2_put_prox <- S2_em[1,]
S2_put_dist <- S2_em[2,]
S2_bg_prox <- S2_em[3,]
S2_bg_dist <- S2_em[4,]


fisher.test(rbind(S2_put_prox, S2_put_dist))
fisher.test(rbind(S2_put_prox, S2_bg_prox))
fisher.test(rbind(S2_put_prox, S2_bg_dist))
fisher.test(rbind(S2_put_dist, S2_bg_prox))
fisher.test(rbind(S2_put_dist, S2_bg_dist))
fisher.test(rbind(S2_bg_prox, S2_bg_dist))

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

pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure3/Figure3wbg.pdf",
width = 28, height = 14, pointsize = 14)
par(mar=c(4,5.5,4,4), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)

  barplot(BG3_per, beside=T, horiz = T, col = cbbPalette[7:4],
      xlim = c(0,100), xlab = "Expressed Contacted Genes", xaxt = "n",
      ylab = "Contact Classification", main = "BG3 Contacted Regions\nWith Expression")
    axis(1, seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"), las=2)
    # putative proximal and distal only
      lines(c(max(BG3_per) + 1.5, max(BG3_per) + 3), c(1.8,1.8))
      lines(c(max(BG3_per) + 3, max(BG3_per) + 3), c(0.8,1.8))
      lines(c(max(BG3_per) + 1.5, max(BG3_per) + 3), c(0.8,0.8))
      text(x = max(BG3_per) + 4, y = 0.8 + (1.8-0.8)/2, "***", srt = 90)
    # putative proximal and bg proximal
      lines(c(max(BG3_per) + 4, max(BG3_per) + 5.5), c(3.0,3.0))
      lines(c(max(BG3_per) + 5.5, max(BG3_per) + 5.5), c(0.8,3.0))
      lines(c(max(BG3_per) + 4, max(BG3_per) + 5.5), c(0.8,0.8))
      text(x = max(BG3_per) + 6.5, y = 0.8 + (3.0-0.8)/2, "***", srt = 90)
    # putative proximal and bg distal only
      lines(c(max(BG3_per) + 6.5, max(BG3_per) + 8), c(4.2,4.2))
      lines(c(max(BG3_per) + 8, max(BG3_per) + 8), c(0.8,4.2))
      lines(c(max(BG3_per) + 6.5, max(BG3_per) + 8), c(0.8,0.8))
      text(x = max(BG3_per) + 9, y = 0.8 + (4.2-0.8)/2, "***", srt = 90)
    # distal only and bg proximal
      lines(c(max(BG3_per) + 1.5, max(BG3_per) + 3), c(3.0,3.0))
      lines(c(max(BG3_per) + 3, max(BG3_per) + 3), c(2.0,3.0))
      lines(c(max(BG3_per) + 1.5, max(BG3_per) + 3), c(2.0,2.0))
      text(x = max(BG3_per) + 4, y = 2.0 + (3.0-2.0)/2, "***", srt = 90)
    # bg proximal and bg distal
      lines(c(max(BG3_per) + 1.5, max(BG3_per) + 3), c(4.2,4.2))
      lines(c(max(BG3_per) + 3, max(BG3_per) + 3), c(3.2,4.2))
      lines(c(max(BG3_per) + 1.5, max(BG3_per) + 3), c(3.2,3.2))
      text(x = max(BG3_per) + 4, y = 3.2 + (4.2-3.2)/2, "***", srt = 90)

  barplot(S2_per, beside=T, horiz = T, col = cbbPalette[7:4],
      xlim = c(0,100), xlab = "Expressed Contacted Genes", xaxt = "n",
      ylab = "Contact Classification", main = "S2 Contacted Regions\nWith Expression")
    axis(1, seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"), las=2)
      # putative proximal and distal only
        lines(c(max(BG3_per) + 1.5, max(BG3_per) + 3), c(1.8,1.8))
        lines(c(max(BG3_per) + 3, max(BG3_per) + 3), c(0.8,1.8))
        lines(c(max(BG3_per) + 1.5, max(BG3_per) + 3), c(0.8,0.8))
        text(x = max(BG3_per) + 4, y = 0.8 + (1.8-0.8)/2, "***", srt = 90)
      # putative proximal and bg proximal
        lines(c(max(BG3_per) + 4, max(BG3_per) + 5.5), c(3.0,3.0))
        lines(c(max(BG3_per) + 5.5, max(BG3_per) + 5.5), c(0.8,3.0))
        lines(c(max(BG3_per) + 4, max(BG3_per) + 5.5), c(0.8,0.8))
        text(x = max(BG3_per) + 6.5, y = 0.8 + (3.0-0.8)/2, "***", srt = 90)
      # putative proximal and bg distal only
        lines(c(max(BG3_per) + 6.5, max(BG3_per) + 8), c(4.2,4.2))
        lines(c(max(BG3_per) + 8, max(BG3_per) + 8), c(0.8,4.2))
        lines(c(max(BG3_per) + 6.5, max(BG3_per) + 8), c(0.8,0.8))
        text(x = max(BG3_per) + 9, y = 0.8 + (4.2-0.8)/2, "***", srt = 90)
      # distal only and bg proximal
        lines(c(max(BG3_per) + 1.5, max(BG3_per) + 3), c(3.0,3.0))
        lines(c(max(BG3_per) + 3, max(BG3_per) + 3), c(2.0,3.0))
        lines(c(max(BG3_per) + 1.5, max(BG3_per) + 3), c(2.0,2.0))
        text(x = max(BG3_per) + 4, y = 2.0 + (3.0-2.0)/2, "***", srt = 90)
      # S2 only distal only and bg distal
        lines(c(max(BG3_per) + 9.5, max(BG3_per) + 11), c(4.2,4.2))
        lines(c(max(BG3_per) + 11, max(BG3_per) + 11), c(2.0,4.2))
        lines(c(max(BG3_per) + 9.5, max(BG3_per) + 11), c(2.0,2.0))
        text(x = max(BG3_per) + 12, y = 2.0 + (4.2-2.0)/2, "***", srt = 90)
      # bg proximal and bg distal
        lines(c(max(BG3_per) + 1.5, max(BG3_per) + 3), c(4.2,4.2))
        lines(c(max(BG3_per) + 3, max(BG3_per) + 3), c(3.2,4.2))
        lines(c(max(BG3_per) + 1.5, max(BG3_per) + 3), c(3.2,3.2))
        text(x = max(BG3_per) + 4, y = 3.2 + (4.2-3.2)/2, "***", srt = 90)

dev.off()
