
# Setting Working Directory
setwd("~/Human_Datasets/H9/csv_files/")

# Colourblind accessible palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Loading in results from csv file
results <- read.csv("Model_Analysis.csv", header = T)

# Converting to a matrix to plot
plot_matrix <- as.matrix(results[,2:9])
rownames(plot_matrix) <- results$X

# Barplots
# Setting colour
cols <- cbbPalette[2:4]

setwd("~/Human_Datasets/H9/Plots/Model_selection_plots")
colnames(plot_matrix) <- c(letters[1:8])

# Average recall stats for whole dataset and testing
pdf(file = "Model_comparisons.pdf",
width = 18, height = 5, pointsize = 14)
par(mar = c(5, 5, 5, 8), xpd = T, cex = 1.2)

# XAI_testing_recalls
barplot(plot_matrix, beside = T, ylim = c(62, 86), col = cols, xpd = F,
main = "Novel Data Recall", yaxt = "none",
xlab = "Model", ylab = "Recall")
axis(2, seq(62, 86, 4), labels = paste0(seq(62, 86, 4), "%"), las=2)
legend(33, 82, c("Average", "Non-\nEnhancer", "Enhancer"), fill = cols)

dev.off()
