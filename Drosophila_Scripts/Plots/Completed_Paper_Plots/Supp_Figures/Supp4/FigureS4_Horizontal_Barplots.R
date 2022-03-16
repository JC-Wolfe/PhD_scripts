
library(GenomicRanges)
library(rtracklayer)

BG3_distal_only_glist <- get(load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/BG3_starr_distal_only_glist.Rda"))
BG3_proximal_enhancers_glist <- get(load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/BG3_starr_proximal_glist.Rda"))

S2_distal_only_glist <- get(load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/S2_starr_distal_only_glist.Rda"))
S2_proximal_enhancers_glist <- get(load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/S2_starr_proximal_glist.Rda"))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Proximal BG3 (BG3P)

# Now which of these don't have any expresion?
no_expression_BG3P <- BG3_proximal_enhancers_glist[sapply(BG3_proximal_enhancers_glist, length) == 0]

# And which of these are expressed?
enhancers_contacted_BG3P <- BG3_proximal_enhancers_glist[sapply(BG3_proximal_enhancers_glist, length) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$BG3_FPKM, na.rm = T))
}

max_to_plot_BG3P <- sapply(enhancers_contacted_BG3P, maxExp)

no_expression_max_BG3P <- c(no_expression_BG3P, enhancers_contacted_BG3P[max_to_plot_BG3P == 0])
expression_max_BG3P <- max_to_plot_BG3P[max_to_plot_BG3P > 0]

# DISTAL ONLY BG3 (BG3D)

# Now which of these don't have any expresion?
no_expression_BG3D <- BG3_distal_only_glist[sapply(BG3_distal_only_glist, length) == 0]

# And which of these are expressed?
enhancers_contacted_BG3D <- BG3_distal_only_glist[sapply(BG3_distal_only_glist, length) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$BG3_FPKM, na.rm = T))
}

max_to_plot_BG3D <- sapply(enhancers_contacted_BG3D, maxExp)

no_expression_max_BG3D <- c(no_expression_BG3D, enhancers_contacted_BG3D[max_to_plot_BG3D == 0])
expression_max_BG3D <- max_to_plot_BG3D[max_to_plot_BG3D > 0]

BG3_em <- matrix(0, 2, 2, dimnames = list(c("Proximal", "Distal Only"),
  c("Expressed", "Not Expressed")))

BG3_em[1,1] <- length(expression_max_BG3P)
BG3_em[1,2] <- length(no_expression_max_BG3P)
BG3_em[2,1] <- length(expression_max_BG3D)
BG3_em[2,2] <- length(no_expression_max_BG3D)

fisher.test(BG3_em)




# Proximal S2 (S2P)

# Now which of these don't have any expresion?
no_expression_S2P <- S2_proximal_enhancers_glist[sapply(S2_proximal_enhancers_glist, length) == 0]

# And which of these are expressed?
enhancers_contacted_S2P <- S2_proximal_enhancers_glist[sapply(S2_proximal_enhancers_glist, length) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$S2_FPKM, na.rm = T))
}

max_to_plot_S2P <- sapply(enhancers_contacted_S2P, maxExp)

no_expression_max_S2P <- c(no_expression_S2P, enhancers_contacted_S2P[max_to_plot_S2P == 0])
expression_max_S2P <- max_to_plot_S2P[max_to_plot_S2P > 0]

# DISTAL ONLY S2 (S2D)

# Now which of these don't have any expresion?
no_expression_S2D <- S2_distal_only_glist[sapply(S2_distal_only_glist, length) == 0]

# And which of these are expressed?
enhancers_contacted_S2D <- S2_distal_only_glist[sapply(S2_distal_only_glist, length) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$S2_FPKM, na.rm = T))
}

max_to_plot_S2D <- sapply(enhancers_contacted_S2D, maxExp)

no_expression_max_S2D <- c(no_expression_S2D, enhancers_contacted_S2D[max_to_plot_S2D == 0])
expression_max_S2D <- max_to_plot_S2D[max_to_plot_S2D > 0]

S2_em <- matrix(0, 2, 2, dimnames = list(c("Proximal", "Distal Only"),
  c("Expressed", "Not Expressed")))

S2_em[1,1] <- length(expression_max_S2P)
S2_em[1,2] <- length(no_expression_max_S2P)
S2_em[2,1] <- length(expression_max_S2D)
S2_em[2,2] <- length(no_expression_max_S2D)

fisher.test(S2_em)

BG3_per <- BG3_em[,1]
BG3_per[1] <- BG3_em[1,1]/(BG3_em[1,1]+BG3_em[1,2])*100
BG3_per[2] <- BG3_em[2,1]/(BG3_em[2,1]+BG3_em[2,2])*100

S2_per <- S2_em[,1]
S2_per[1] <- S2_em[1,1]/(S2_em[1,1]+S2_em[1,2])*100
S2_per[2] <- S2_em[2,1]/(S2_em[2,1]+S2_em[2,2])*100

pdf(file = "/home/jw18713/Project1/Paper_Plots/FigureS4/FigureS4C.pdf",
width = 14, height = 5, pointsize = 14)
par(mar=c(4,5.5,4,4), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)

  barplot(BG3_per, beside=T, horiz = T, col = cbbPalette[7:6],
      xlim = c(0,100), xlab = "Expressed Contacted Genes\n(STARR-seq Only Enhancers)", xaxt = "n",
      ylab = "Contact Classification", main = "BG3 Contacted Regions\nWith Expression")
    axis(1, seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"), las=2)
    lines(c(max(BG3_per) + 2.5, max(BG3_per) + 5), c(1.9,1.9))
    lines(c(max(BG3_per) + 5, max(BG3_per) + 5), c(0.7,1.9))
    lines(c(max(BG3_per) + 2.5, max(BG3_per) + 5), c(0.7,0.7))
    text(x = max(BG3_per) + 7.5, y = 0.7 + (1.9-0.7)/2, "***", srt = 90)

  barplot(S2_per, beside=T, horiz = T, col = cbbPalette[7:6],
      xlim = c(0,100), xlab = "Expressed Contacted Genes\n(STARR-seq Only Enhancers)", xaxt = "n",
      ylab = "Contact Classification", main = "S2 Contacted Regions\nWith Expression")
    axis(1, seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"), las=2)
    lines(c(max(S2_per) + 2.5, max(S2_per) + 5), c(1.9,1.9))
    lines(c(max(S2_per) + 5, max(S2_per) + 5), c(0.7,1.9))
    lines(c(max(S2_per) + 2.5, max(S2_per) + 5), c(0.7,0.7))
    text(x = max(S2_per) + 7.5, y = 0.7 + (1.9-0.7)/2, "***", srt = 90)

dev.off()
