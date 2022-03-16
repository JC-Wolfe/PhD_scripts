
# Setting Working Directory
setwd("Human_Datasets/H9/csv_files/")

# Colourblind accessible palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7")

# Loading in results from csv file
results <- read.csv("25%_rule_number_recall.csv")

# Matrix creation for average recall comparison plot
WD_averages <- matrix(results$WD_average, nrow = 1, ncol = 5,
  dimnames = list(c("Recall"), c("450", "250", "100", "50", "25")))

Test_averages <- matrix(results$Testing_average, nrow = 1, ncol = 5,
  dimnames = list(c("Recall"), c("450", "250", "100", "50", "25")))

# Matrix creation for balance plots
WD_balance <- matrix(c(results$WD_non_enhancer, results$WD_enhancer),
  nrow = 2, ncol = 5, dimnames = list(c("Non-Enhancer", "Enhancer"),
  c("450", "250", "100", "50", "25")), byrow = T)

Test_balance <- matrix(c(results$Testing_non_enhancer, results$Testing_enhancer),
  nrow = 2, ncol = 5, dimnames = list(c("Non-Enhancer", "Enhancer"),
  c("450", "250", "100", "50", "25")), byrow = T)

# Barplots
# Setting colourl
cols <- cbbPalette[2:3]

setwd("~/Human_Datasets/H9/Plots/Model_selection_plots")

# Average recall stats for whole dataset and testing
pdf(file = "Rule_number_recall_stats.pdf",
width = 18, height = 11, pointsize = 14)
par(mar = c(5, 5, 5, 9), xpd = T, cex = 1.2, mfrow = c(2,2))

# Testing_averages
barplot(Test_averages, beside = T, ylim = c(70, 80), col = cols[2], xpd = F,
main = "Testing Average Recall", yaxt = "none",
xlab = "Number of Rules", ylab = "Average Recall")
axis(2, seq(70, 80, 2), labels = paste0(seq(70, 80, 2), "%"), las=2)

# WD_averages
barplot(WD_averages, beside = T, ylim = c(70, 80), col = cols[2], xpd = F,
main = "Whole Dataset Average Recall", yaxt = "none",
xlab = "Number of Rules", ylab = "Averages Recall")
axis(2, seq(70, 80, 2), labels = paste0(seq(70, 80, 2), "%"), las=2)

# testing_balance
barplot(Test_balance, beside = T, ylim = c(65, 85), col = cols, xpd = F,
main = "Testing Recall Balance", yaxt = "none",
xlab = "Number of Rules", ylab = "Recall")
axis(2, seq(65, 85, 5), labels = paste0(seq(65, 85, 5), "%"), las=2)
legend(16, 72.5, c("Non-\nEnhancer", "Enhancer"), fill = cols)

# WD_balance
barplot(WD_balance, beside = T, ylim = c(65, 85), col = cols, xpd = F,
main = "Whole Dataset Recall Balance", yaxt = "none",
xlab = "Number of Rules", ylab = "Recall")
axis(2, seq(65, 85, 5), labels = paste0(seq(65, 85, 5), "%"), las=2)
legend(16, 72.5, c("Non-\nEnhancer", "Enhancer"), fill = cols)

dev.off()
