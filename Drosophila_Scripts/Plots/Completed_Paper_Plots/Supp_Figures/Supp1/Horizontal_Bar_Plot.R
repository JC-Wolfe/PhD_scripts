
pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure3/Figure3_Horiz_BG.pdf",
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
    # putative proximal and proximal bg
      lines(c(max(BG3_plot[,2]) + 1.5, max(BG3_plot[,2]) + 3), c(2.5,2.5))
      lines(c(max(BG3_plot[,2]) + 3, max(BG3_plot[,2]) + 3), c(1.5,2.5))
      lines(c(max(BG3_plot[,2]) + 1.5, max(BG3_plot[,2]) + 3), c(1.5,1.5))
      text(x = max(BG3_plot[,2]) + 4, y = 1.5 + (2.5-1.5)/2, "***", srt = 90)
    # putative distal only and putative proximal
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
    # putative proximal and proximal bg
      lines(c(max(S2_plot[,2]) + 1.5, max(S2_plot[,2]) + 3), c(2.5,2.5))
      lines(c(max(S2_plot[,2]) + 3, max(S2_plot[,2]) + 3), c(1.5,2.5))
      lines(c(max(S2_plot[,2]) + 1.5, max(S2_plot[,2]) + 3), c(1.5,1.5))
      text(x = max(S2_plot[,2]) + 4, y = 1.5 + (2.5-1.5)/2, "**", srt = 90)
    # putative distal only and distal bg
      lines(c(max(S2_plot[,1]) + 1.5, max(S2_plot[,1]) + 3), c(4.5,4.5))
      lines(c(max(S2_plot[,1]) + 3, max(S2_plot[,1]) + 3), c(5.5,4.5))
      lines(c(max(S2_plot[,1]) + 1.5, max(S2_plot[,1]) + 3), c(5.5,5.5))
      text(x = max(S2_plot[,1]) + 4, y = 4.5 + (5.5-4.5)/2, "***", srt = 90)
    # putative distal only and putative proximal
      lines(c(max(S2_plot[,1]) + 4.5, max(S2_plot[,1]) + 6), c(5.5,5.5))
      lines(c(max(S2_plot[,1]) + 6, max(S2_plot[,1]) + 6), c(2.5,5.5))
      lines(c(max(S2_plot[,1]) + 4.5, max(S2_plot[,1]) + 6), c(2.5,2.5))
      text(x = max(S2_plot[,1]) + 7, y = 2.5 + (5.5-2.5)/2, "***", srt = 90)

dev.off()
