library(GenomicRanges)
library(rtracklayer)

load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/BG3_distal_only_glist.Rda")
load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/BG3_proximal_enhancers_glist.Rda")

load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/S2_distal_only_glist.Rda")
load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/S2_proximal_enhancer_glist.Rda")

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

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Ok, so this is my histogram then
pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure3/Figure3I&J.pdf",
width = 14, height = 6, pointsize = 14)
par(mar=c(4.9,5.5,5.5,8), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)

hist(log2(expression_max_BG3D+1), main="BG3 Enhancers\nMax Contacted Expression",
  xlab="Expression FPKM", col="#0072B280", cex.main=2, cex.lab=1.5, axes=F,
  ylim = c(0,2500), breaks = 16)
axis(1, at=seq(0,16,by=4), labels=c(0,2^seq(2,16, by = 4)))
axis(2, at=seq(0,2500, by = 500), labels = seq(0,2500, by = 500))

hist(log2(expression_max_BG3P+1), col="#D55E0080", add=T, breaks = 16)

legend(16.5,1500,c("Distal", "Proximal"), fill=c("#0072B280", "#D55E0080"))

# S2 stuff

hist(log2(expression_max_S2D+1), main="S2 Enhancers\nMax Contacted Expression",
  xlab="Expression FPKM", col="#0072B280", cex.main=2, cex.lab=1.5, axes=F,
  ylim = c(0,2500), breaks = 16)
axis(1, at=seq(0,16,by=4), labels=c(0,2^seq(4,16, by = 4)))
axis(2, at=seq(0,2500, by = 500), labels = seq(0,2500, by = 500))

hist(log2(expression_max_S2P+1), col="#D55E0080", add=T, breaks = 16)

legend(16.5,1500,c("Distal", "Proximal"), fill=c("#0072B280", "#D55E0080"))

dev.off()
