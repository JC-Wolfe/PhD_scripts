
setwd("~/Pre_viva/Hi-C/grange_lists")

load("H9_background_distal_EA.Rda")
load("H9_background_proximal_EA.Rda")
load("H9_common_distal_EA.Rda")
load("H9_common_proximal_EA.Rda")
load("H9_putative_distal_EA.Rda")
load("H9_putative_proximal_EA.Rda")
load("H9_starr_only_distal_EA.Rda")
load("H9_starr_only_proximal_EA.Rda")
load("Neither_matrix.Rda")

# Calculating all with contacts over all total (with and without)

all_calc <- function(x, y, neither){
  # Calculating all with contacts over all total (with and without)
  all <- sum(x$count_with_contacts, y$count_with_contacts)/
    sum(x$count_without_contacts, y$count_without_contacts,
    x$count_with_contacts,y$count_with_contacts, neither)
  return(all)
}

H9_matrix <- matrix(0, 4, 3, dimnames = list(c("Putative", "Common",
  "STARR-seq only", "Whole Genome"), c("Proximal", "Distal Only", "All")))

H9_matrix[1,1] <- H9_putative_proximal_EA$contact_proportion
H9_matrix[1,2] <- H9_putative_proximal_EA$contact_proportion
H9_matrix[1,3] <- all_calc(H9_putative_proximal_EA, H9_putative_proximal_EA,
  neither = Neither_matrix[1,1])

H9_matrix[2,1] <- H9_common_proximal_EA$contact_proportion
H9_matrix[2,2] <- H9_common_distal_EA$contact_proportion
H9_matrix[2,3] <- all_calc(H9_common_proximal_EA, H9_common_distal_EA,
  neither = Neither_matrix[1,2])

H9_matrix[3,1] <- H9_starr_only_proximal_EA$contact_proportion
H9_matrix[3,2] <- H9_starr_only_distal_EA$contact_proportion
H9_matrix[3,3] <- all_calc(H9_starr_only_proximal_EA, H9_starr_only_distal_EA,
  neither = Neither_matrix[1,3])

H9_matrix[4,1] <- H9_background_proximal_EA$contact_proportion
H9_matrix[4,2] <- H9_background_distal_EA$contact_proportion
H9_matrix[4,3] <- all_calc(H9_background_proximal_EA, H9_background_distal_EA,
  neither = Neither_matrix[1,4])


setwd("~/Pre_viva/Hi-C/Plots")

#---------------------------------H9 Stats (All)------------------------------#
#------------------------------------------------------------------------------#
all_counts <- function(x, y, neither){
  without <- sum(x$count_without_contacts, y$count_without_contacts, neither)
  with <- sum(x$count_with_contacts, y$count_with_contacts)
  return(c(without = without, with = with))
}

# Putative
H9_putative_all <- all_counts(H9_putative_proximal_EA, H9_putative_proximal_EA,
  neither = Neither_matrix[1,1])
# Common
H9_common_all <- all_counts(H9_common_proximal_EA, H9_common_distal_EA,
  neither = Neither_matrix[1,2])
# STARR-seq Only
H9_starr_all <- all_counts(H9_starr_only_proximal_EA, H9_starr_only_distal_EA,
  neither = Neither_matrix[1,3])
# Background
H9_background_all <- all_counts(H9_background_proximal_EA, H9_background_distal_EA,
  neither = Neither_matrix[1,4])

H9_all_stats <- matrix("", 4, 4, dimnames = list(c(rownames(H9_matrix)),
rownames(H9_matrix)))

same <- matrix(c(1,1,2,2,3,3,4,4), 4, 2 , byrow = T)

# Filling putative/putative etc with -
H9_all_stats[same] <- "-"

# Putative/Common
H9_all_stats[2,1] <- fisher.test(matrix(c(H9_putative_all, H9_common_all),
  ncol = 2, byrow = T))$p.value

# Putative/STARR-seq only
H9_all_stats[3,1] <- fisher.test(matrix(c(H9_putative_all, H9_starr_all),
  ncol = 2, byrow = T))$p.value

# Putative/Whole Genome
H9_all_stats[4,1] <- fisher.test(matrix(c(H9_putative_all, H9_background_all),
  ncol = 2, byrow = T))$p.value

# Common/STARR-seq only
H9_all_stats[3,2] <- fisher.test(matrix(c(H9_common_all, H9_starr_all),
  ncol = 2, byrow = T))$p.value

# Common/Whole Genome
H9_all_stats[4,2] <- fisher.test(matrix(c(H9_common_all, H9_background_all),
  ncol = 2, byrow = T))$p.value

# STARR-seq/Whole Genome
H9_all_stats[4,3] <- fisher.test(matrix(c(H9_starr_all, H9_background_all),
  ncol = 2, byrow = T))$p.value

# Writing stats matrix to a csv file
write.csv(H9_all_stats,
  file = "Stats/Figure3C_All.csv")


#---------------------------H9 Stats (Distal Only)----------------------------#
#------------------------------------------------------------------------------#

# Putative
H9_putative_d <- c(H9_putative_proximal_EA$count_without_contacts,
  H9_putative_proximal_EA$count_with_contacts)
# Common
H9_common_d <- c(H9_common_distal_EA$count_without_contacts,
  H9_common_distal_EA$count_with_contacts)
# STARR-seq Only
H9_starr_d <- c(H9_starr_only_distal_EA$count_without_contacts,
  H9_starr_only_distal_EA$count_with_contacts)
# Background
H9_background_d <- c(H9_background_distal_EA$count_without_contacts,
  H9_background_distal_EA$count_with_contacts)

H9_d_stats <- matrix("", 4, 4, dimnames = list(c(rownames(H9_matrix)),
rownames(H9_matrix)))

same <- matrix(c(1,1,2,2,3,3,4,4), 4, 2 , byrow = T)

# Filling putative/putative etc with -
H9_d_stats[same] <- "-"

# Putative/Common
H9_d_stats[2,1] <- fisher.test(matrix(c(H9_putative_d, H9_common_d),
  ncol = 2, byrow = T))$p.value

# Putative/STARR-seq only
H9_d_stats[3,1] <- fisher.test(matrix(c(H9_putative_d, H9_starr_d),
  ncol = 2, byrow = T))$p.value

# Putative/Whole Genome
H9_d_stats[4,1] <- fisher.test(matrix(c(H9_putative_d, H9_background_d),
  ncol = 2, byrow = T))$p.value

# Common/STARR-seq only
H9_d_stats[3,2] <- fisher.test(matrix(c(H9_common_d, H9_starr_d),
  ncol = 2, byrow = T))$p.value

# Common/Whole Genome
H9_d_stats[4,2] <- fisher.test(matrix(c(H9_common_d, H9_background_d),
  ncol = 2, byrow = T))$p.value

# STARR-seq/Whole Genome
H9_d_stats[4,3] <- fisher.test(matrix(c(H9_starr_d, H9_background_d),
  ncol = 2, byrow = T))$p.value

# Writing stats matrix to a csv file
write.csv(H9_d_stats,
  file = "Stats/Figure3C_Distal_Only.csv")


#-----------------------------H9 Stats (Proximal)-----------------------------#
#------------------------------------------------------------------------------#

# Putative
H9_putative_pr <- c(H9_putative_proximal_EA$count_without_contacts,
  H9_putative_proximal_EA$count_with_contacts)
# Common
H9_common_pr <- c(H9_common_proximal_EA$count_without_contacts,
  H9_common_proximal_EA$count_with_contacts)
# STARR-seq Only
H9_starr_pr <- c(H9_starr_only_proximal_EA$count_without_contacts,
  H9_starr_only_proximal_EA$count_with_contacts)
# Background
H9_background_pr <- c(H9_background_proximal_EA$count_without_contacts,
  H9_background_proximal_EA$count_with_contacts)

H9_pr_stats <- matrix("", 4, 4, dimnames = list(c(rownames(H9_matrix)),
rownames(H9_matrix)))

same <- matrix(c(1,1,2,2,3,3,4,4), 4, 2 , byrow = T)

# Filling putative/putative etc with -
H9_pr_stats[same] <- "-"

# Putative/Common
H9_pr_stats[2,1] <- fisher.test(matrix(c(H9_putative_pr, H9_common_pr),
  ncol = 2, byrow = T))$p.value

# Putative/STARR-seq only
H9_pr_stats[3,1] <- fisher.test(matrix(c(H9_putative_pr, H9_starr_pr),
  ncol = 2, byrow = T))$p.value

# Putative/Whole Genome
H9_pr_stats[4,1] <- fisher.test(matrix(c(H9_putative_pr, H9_background_pr),
  ncol = 2, byrow = T))$p.value

# Common/STARR-seq only
H9_pr_stats[3,2] <- fisher.test(matrix(c(H9_common_pr, H9_starr_pr),
  ncol = 2, byrow = T))$p.value

# Common/Whole Genome
H9_pr_stats[4,2] <- fisher.test(matrix(c(H9_common_pr, H9_background_pr),
  ncol = 2, byrow = T))$p.value

# STARR-seq/Whole Genome
H9_pr_stats[4,3] <- fisher.test(matrix(c(H9_starr_pr, H9_background_pr),
  ncol = 2, byrow = T))$p.value

# Writing stats matrix to a csv file
write.csv(H9_pr_stats,
  file = "Stats/Figure3C_Proximal.csv")


#-------------------------------------Plotting---------------------------------#
#------------------------------------------------------------------------------#

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
  "#D55E00", "#CC79A7")
cols <- c("white", cbPalette[c(5,2,7)])


pdf(file = "Figure3C.pdf", width = 10,
  height = 6, pointsize = 14)
par(mar=c(4,3.5,4,13), cex = 1.2)
par(xpd = NA)

  barplot(H9_matrix[4:1,], beside=T, horiz = T, col = cols,
      xlim = c(0,1), xlab = "% of enhancers that contact expressed genes", xaxt = "n",
      main = "Enhancers in H9 contact expressed genes")
    axis(1, seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10), "%"), las=2)
    legend(1.1,12.5,rownames(H9_matrix), fill=rev(cols), bty='n')

dev.off()
