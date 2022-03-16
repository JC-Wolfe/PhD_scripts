# Importing libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6")
BiocManager::install("genomation")

library(BSgenome.Dmelanogaster.UCSC.dm6)
library(GenomicRanges)
library(rtracklayer)
library(genomation)

setwd("~/P1_Full_Run/Datasets/S2/S2_STARR-seq/STARK2")

# STARK datasets
GEO_txt_convert <- function(x){
  c_tab <- read.table(x)
  c_gr <- GRanges(seqnames = c_tab[,1],
    ranges = IRanges(start = c_tab[,2] - 200, end = c_tab[,2] + 200))
  return(c_gr)
}

# Developmental
d1 <- GEO_txt_convert("GSE57876_dCP_S2_merged.peaks.txt")
d2 <- GEO_txt_convert("GSE57876_eEF1delta_S2_rep1.peaks.txt")
dCP <- reduce(c(d1, d2))

# Housekeeping
hk1 <- GEO_txt_convert("GSE57876_hkCP_S2_merged.peaks.txt")
hk2 <- GEO_txt_convert("GSE57876_NipB_S2_rep2.peaks.txt")
hk3 <- GEO_txt_convert("GSE57876_pnr_S2.peaks.txt")
hk4 <- GEO_txt_convert("GSE57876_X16_S2_rep2.peaks.txt")
hkCP <- reduce(c(hk1, hk2, hk3, hk4))

# Getting the dm6 genome
genome <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6)

# Choosing a bin size
binsize <- 10

# Creation of the GRange object to contain ChIP metadata
chrs <- seqnames(Dmelanogaster)[1:5]
S2_tiled <- GRanges()
seqlevels(S2_tiled) <- chrs
for(i in 1:length(chrs)){
  print(chrs[i])
  start <- seq(1, length(Dmelanogaster[[i]]), by=binsize)
  buffer <- GRanges(seqnames=chrs[i], ranges=IRanges(start = start, end = start + 9), strand="*")
  S2_tiled <- c(S2_tiled, buffer)
}

# Enhancers to cover
d_only <- setdiff(dCP, hkCP)
hk_only <- setdiff(hkCP, dCP)
coverage <- c(d_only, hk_only)
overlaps <- findOverlaps(S2_tiled, coverage)
enhancers_tiled <- S2_tiled[queryHits(overlaps)]


#Setting the right working directory
setwd("/home/jw18713/P1_Full_Run/Datasets/S2/S2_histone")

# Creating a list of all of the folders that I need (excluding non .wig files)
folders <- dir()

# Using paste0 to create a list of file paths to folders containing each .wig file
paths <- paste0(getwd(),"/",folders,"/signal_data_files/")

# Creating a list of "0"s to be filled with the relevant files and file names
files <- rep("0", length(folders))
name <- rep("0", length(folders))

# A for loop to loop over every directory contained in paths that:
# Sets the working directory to the current element in paths
# Locates the smoothed file using grep and saves it to a temporary variable
# Changes the "i"th element of names to the correct dataset name
# Changes the "i"th element of files to the path to the relevant smoothed .wig file
for(i in seq_along(paths)){
setwd(paths[i])
smooth <- dir()[grep("smoothedM",dir())]
name[i] <- strsplit(smooth,":")[[1]][1]
files[i] <- paste0(paths[i],smooth)
}


# Chain for lift over
chain <- import.chain("/home/jw18713/P1_Full_Run/Misc_Resources/dm3ToDm6.over.chain")

for(i in seq_along(files)){
  gro <- readGeneric(files[i], skip=1, sep=" ", meta.cols=list(score=4))
  seqlevelsStyle(gro) <- "UCSC"
  gro <- unlist(liftOver(gro, chain))
  # WE HAVE LIFT OVER!
  overlaps <- findOverlaps(enhancers_tiled,gro)
  scores <- rep(0, length(enhancers_tiled))
  scores[as.vector(queryHits(overlaps))] <- gro$score[as.vector(subjectHits(overlaps))]
  print(name[i])
  mcols(enhancers_tiled)[i] <- scores
  names(mcols(enhancers_tiled))[i] <- name[i]
}


sortdf <- mcols(enhancers_tiled)[,sort(names(mcols(enhancers_tiled)))]
mcols(enhancers_tiled) <- sortdf

# Enhancer labelling (Developmental minority class)
scoring <- rep(0, length(enhancers_tiled))
overlaps <- findOverlaps(d_only, enhancers_tiled)
scoring[subjectHits(overlaps)] <- 1
enhancers_tiled$Developmental <- scoring

# Getting min and max values to de-normalise outputs later
output_max <- max(mcols(enhancers_tiled)[length(mcols(enhancers_tiled))][,1])
output_min <- min(mcols(enhancers_tiled)[length(mcols(enhancers_tiled))][,1])
save(output_max, file="/home/jw18713/P1_Full_Run/Rda_objects/Denorm_data/S2_dCP_hkCP_dC_denorm_max.Rda")
save(output_min, file="/home/jw18713/P1_Full_Run/Rda_objects/Denorm_data/S2_dCP_hkCP_denorm_min.Rda")

# Scaling the dataset
for (c in seq_along(mcols(enhancers_tiled))){
  max_value = max(mcols(enhancers_tiled)[c][,1])
  min_value = min(mcols(enhancers_tiled)[c][,1])
  mcols(enhancers_tiled)[c][,1] <- (mcols(enhancers_tiled)[c][,1] - min_value) / (max_value - min_value)
}


save(enhancers_tiled, file="/home/jw18713/P1_Full_Run/Rda_objects/Tiled_genomes/S2_dCP_hkCP_tiled.Rda")

dfwrite <- as.data.frame(enhancers_tiled)
write.csv(dfwrite, file="/home/jw18713/P1_Full_Run/csv_files/S2_dCP_hkCP_tiled.csv", row.names=F)
