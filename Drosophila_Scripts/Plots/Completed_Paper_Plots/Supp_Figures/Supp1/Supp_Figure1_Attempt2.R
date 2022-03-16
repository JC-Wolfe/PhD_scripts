install.packages("PRROC")
library(PRROC)

BG3_XAI <- get(load("/home/jw18713/Project1/data/model_predictions/BG3_FuzzyGR.Rda"))
S2_XAI <- get(load("/home/jw18713/Project1/data/model_predictions/S2_FuzzyGR.Rda"))

fg1 <- BG3_XAI$Conf..Perc.1[BG3_XAI$STARR_seq_binary == 1]
bg1 <- BG3_XAI$Conf..Perc.1[BG3_XAI$STARR_seq_binary == 0]

# BG3 ROC Curve
BG3_roc <- roc.curve(scores.class0 = fg1, scores.class1 = bg1, curve = T)

# BG3 PR Curve
BG3_pr <- pr.curve(scores.class0 = fg1, scores.class1 = bg1, curve = T)

fg2 <- S2_XAI$Conf..Perc.1[S2_XAI$STARR_seq_binary == 1]
bg2 <- S2_XAI$Conf..Perc.1[S2_XAI$STARR_seq_binary == 0]

# S2 ROC Curve
S2_roc <- roc.curve(scores.class0 = fg2, scores.class1 = bg2, curve = T)

# S2 PR Curve
S2_pr <- pr.curve(scores.class0 = fg2, scores.class1 = bg2, curve = T)






pdf(file = "/home/jw18713/Project1/Paper_Plots/FigureS1/Supp_Figure1All.pdf",
width = 12, height = 12, pointsize = 14)
par(mar=c(4,5.5,4,6), cex = 1.2)



plot(BG3_roc, main = "BG3 ROC")
plot(BG3_pr, main = "BG3 PRC")
plot(S2_roc, main = "S2 ROC")
plot(S2_pr, main = "S2 PRC")

dev.off()
