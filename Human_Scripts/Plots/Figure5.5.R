
load("~/Human_Datasets/H1/Processed_datasets/H1_reduction_matrix.Rda")
load("~/Human_Datasets/H9/Processed_datasets/H9_reduction_matrix.Rda")

pdf("~/Human_Datasets/Plots/Supp/Reduction_counts.pdf")
plot(1:7, H9_reduction_matrix[1,], type = "o", main = "Enhancer Counts by\nMinimum Gap Width",
  xlab = "Gap Width", ylab = "Enhancer Count", lty = 1, axes = F,
  col = "blue", ylim = c(0, 3500000))
  points(1:7, H1_reduction_matrix[1,], col = "red", pch = "+")
  lines(1:7, H1_reduction_matrix[1,], col = "red", lty = 2)
  axis(1, at=1:7, labels = colnames(H9_reduction_matrix))
  axis(2, at=seq(0,3500000,by=500000), labels = c("0", "500k", "1M", "1.5M", "2M", "2.5M", "3M", "3.5M"), las = 1)
  legend("topright", legend = c("H9", "H1"), col = c("blue", "red"),
  pch = c("o", "+"), bty = "n")
dev.off()

pdf("~/Human_Datasets/Plots/Supp/Reduction_averages.pdf")
plot(1:7, H9_reduction_matrix[2,], type = "o", main = "Confidence Averages by\nMinimum Gap Width",
  xlab = "Gap Width", ylab = "Confidence Average", lty = 1, axes = F,
  col = "blue", ylim = c(0, 1))
  points(1:7, H1_reduction_matrix[2,], col = "red", pch = "+")
  lines(1:7, H1_reduction_matrix[2,], col = "red", lty = 2)
  axis(1, at=1:7, labels = colnames(H9_reduction_matrix))
  axis(2, at=seq(0,1,by=0.2), labels = seq(0,1,by=0.2), las = 1)
  legend("bottomleft", legend = c("H9", "H1"), col = c("blue", "red"),
  pch = c("o", "+"), bty = "n")
dev.off()

pdf("~/Human_Datasets/Plots/Supp/Reduction_widths.pdf")
plot(1:7, H9_reduction_matrix[3,], type = "o", main = "Enhancer Widths by\nMinimum Gap Width",
  xlab = "Gap Width", ylab = "Average Width", lty = 1, axes = F,
  col = "blue", ylim = c(0, 6000))
  points(1:7, H1_reduction_matrix[3,], col = "red", pch = "+")
  lines(1:7, H1_reduction_matrix[3,], col = "red", lty = 2)
  axis(1, at=1:7, labels = colnames(H9_reduction_matrix))
  axis(2, at=seq(0,6000,by=1500), labels = seq(0,6000,by=1500), las = 1)
  legend("topleft", legend = c("H9", "H1"), col = c("blue", "red"),
  pch = c("o", "+"), bty = "n")
dev.off()
