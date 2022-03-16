library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Loading in BG3 Datasets
BG3_XAi <- get(load("/home/jw18713/Project1/data/model_predictions/BG3_FuzzyGR.Rda"))
BG3_NN <- get(load("/home/jw18713/Project1/data/model_predictions/BG3_NN_GR.Rda"))
BG3_LR <- get(load("/home/jw18713/Project1/data/model_predictions/BG3_LR_GR.Rda"))

# Loading in S2 Datasets
S2_XAi <- get(load("/home/jw18713/Project1/data/model_predictions/S2_FuzzyGR.Rda"))
S2_NN <- get(load("/home/jw18713/Project1/data/model_predictions/S2_NN_GR.Rda"))
S2_LR <- get(load("/home/jw18713/Project1/data/model_predictions/S2_LR_GR.Rda"))

# Classifications for BG3 Fuzzy Logic
TP_BG3_XAi <- length(BG3_XAi[BG3_XAi$STARR_seq_binary == 1 & BG3_XAi$Predicted.Class == 1])
TN_BG3_XAi <- length(BG3_XAi[BG3_XAi$STARR_seq_binary == 0 & BG3_XAi$Predicted.Class == 0])
FP_BG3_XAi <- length(BG3_XAi[BG3_XAi$STARR_seq_binary == 0 & BG3_XAi$Predicted.Class == 1])
FN_BG3_XAi <- length(BG3_XAi[BG3_XAi$STARR_seq_binary == 1 & BG3_XAi$Predicted.Class == 0])

# Classifications for BG3 Neural Network
TP_BG3_NN <- length(BG3_NN[BG3_NN$STARR_seq_binary == 1 & BG3_NN$Predicted.Class == 1])
TN_BG3_NN <- length(BG3_NN[BG3_NN$STARR_seq_binary == 0 & BG3_NN$Predicted.Class == 0])
FP_BG3_NN <- length(BG3_NN[BG3_NN$STARR_seq_binary == 0 & BG3_NN$Predicted.Class == 1])
FN_BG3_NN <- length(BG3_NN[BG3_NN$STARR_seq_binary == 1 & BG3_NN$Predicted.Class == 0])

# Classifications for BG3 Logistic Regression
TP_BG3_LR <- length(BG3_LR[BG3_LR$STARR_seq_binary == 1 & BG3_LR$Predicted.Class == 1])
TN_BG3_LR <- length(BG3_LR[BG3_LR$STARR_seq_binary == 0 & BG3_LR$Predicted.Class == 0])
FP_BG3_LR <- length(BG3_LR[BG3_LR$STARR_seq_binary == 0 & BG3_LR$Predicted.Class == 1])
FN_BG3_LR <- length(BG3_LR[BG3_LR$STARR_seq_binary == 1 & BG3_LR$Predicted.Class == 0])

# Classifications for S2
TP_S2_XAi <- length(S2_XAi[S2_XAi$STARR_seq_binary == 1 & S2_XAi$Predicted.Class == 1])
TN_S2_XAi <- length(S2_XAi[S2_XAi$STARR_seq_binary == 0 & S2_XAi$Predicted.Class == 0])
FP_S2_XAi <- length(S2_XAi[S2_XAi$STARR_seq_binary == 0 & S2_XAi$Predicted.Class == 1])
FN_S2_XAi <- length(S2_XAi[S2_XAi$STARR_seq_binary == 1 & S2_XAi$Predicted.Class == 0])

TP_S2_NN <- length(S2_NN[S2_NN$STARR_seq_binary == 1 & S2_NN$Predicted.Class == 1])
TN_S2_NN <- length(S2_NN[S2_NN$STARR_seq_binary == 0 & S2_NN$Predicted.Class == 0])
FP_S2_NN <- length(S2_NN[S2_NN$STARR_seq_binary == 0 & S2_NN$Predicted.Class == 1])
FN_S2_NN <- length(S2_NN[S2_NN$STARR_seq_binary == 1 & S2_NN$Predicted.Class == 0])

TP_S2_LR <- length(S2_LR[S2_LR$STARR_seq_binary == 1 & S2_LR$Predicted.Class == 1])
TN_S2_LR <- length(S2_LR[S2_LR$STARR_seq_binary == 0 & S2_LR$Predicted.Class == 0])
FP_S2_LR <- length(S2_LR[S2_LR$STARR_seq_binary == 0 & S2_LR$Predicted.Class == 1])
FN_S2_LR <- length(S2_LR[S2_LR$STARR_seq_binary == 1 & S2_LR$Predicted.Class == 0])

# BG3 matrix
BG3_cmat <- matrix(0,3,3)
colnames(BG3_cmat) <- c("Accuracy", "Precision", "Recall")
rownames(BG3_cmat) <- c("XAI", "NN", "LR")

BG3_cmat[1,1] <- TP_BG3_XAi/(TP_BG3_XAi + FN_BG3_XAi)
BG3_cmat[1,2] <- TP_BG3_XAi/(TP_BG3_XAi + FP_BG3_XAi)
BG3_cmat[1,3] <- (TP_BG3_XAi + TN_BG3_XAi)/(TP_BG3_XAi + FP_BG3_XAi + TN_BG3_XAi + FN_BG3_XAi)

BG3_cmat[2,1] <- TP_BG3_NN/(TP_BG3_NN + FN_BG3_NN)
BG3_cmat[2,2] <- TP_BG3_NN/(TP_BG3_NN + FP_BG3_NN)
BG3_cmat[2,3] <- (TP_BG3_NN + TN_BG3_NN)/(TP_BG3_NN + FP_BG3_NN + TN_BG3_NN + FN_BG3_NN)

BG3_cmat[3,1] <- TP_BG3_LR/(TP_BG3_LR + FN_BG3_LR)
BG3_cmat[3,2] <- TP_BG3_LR/(TP_BG3_LR + FP_BG3_LR)
BG3_cmat[3,3] <- (TP_BG3_LR + TN_BG3_LR)/(TP_BG3_LR + FP_BG3_LR + TN_BG3_LR + FN_BG3_LR)

# S2 Matrix

S2_cmat <- matrix(0,3,3)
colnames(S2_cmat) <- c("Accuracy", "Precision", "Recall")
rownames(S2_cmat) <- c("XAI", "NN", "LR")

S2_cmat[1,1] <- TP_S2_XAi/(TP_S2_XAi + FN_S2_XAi)
S2_cmat[1,2] <- TP_S2_XAi/(TP_S2_XAi + FP_S2_XAi)
S2_cmat[1,3] <- (TP_S2_XAi + TN_S2_XAi)/(TP_S2_XAi + FP_S2_XAi + TN_S2_XAi + FN_S2_XAi)

S2_cmat[2,1] <- TP_S2_NN/(TP_S2_NN + FN_S2_NN)
S2_cmat[2,2] <- TP_S2_NN/(TP_S2_NN + FP_S2_NN)
S2_cmat[2,3] <- (TP_S2_NN + TN_S2_NN)/(TP_S2_NN + FP_S2_NN + TN_S2_NN + FN_S2_NN)

S2_cmat[3,1] <- TP_S2_LR/(TP_S2_LR + FN_S2_LR)
S2_cmat[3,2] <- TP_S2_LR/(TP_S2_LR + FP_S2_LR)
S2_cmat[3,3] <- (TP_S2_LR + TN_S2_LR)/(TP_S2_LR + FP_S2_LR + TN_S2_LR + FN_S2_LR)

# Plotting BG3
pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure1/Figure1B&C.pdf",
width = 10, height = 5, pointsize = 14)
par(mar=c(4,5.5,4,6), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)
    barplot(BG3_cmat, col = cbbPalette[2:4], beside = T, ylim = c(0,1),
            main = "BG3 Confusion Matrix Statistics"
            )
    barplot(S2_cmat, col = cbbPalette[2:4], beside = T, ylim = c(0,1),
            main = "S2 Confusion Matrix Statistics"
            )
            legend(13.5,0.75,rownames(BG3_cmat), fill=cbbPalette[2:4])
dev.off()
