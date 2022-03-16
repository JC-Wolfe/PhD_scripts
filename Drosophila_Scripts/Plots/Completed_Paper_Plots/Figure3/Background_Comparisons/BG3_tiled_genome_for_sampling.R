# Importing libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6")
BiocManager::install("genomation")

library(BSgenome.Dmelanogaster.UCSC.dm6)
library(GenomicRanges)
library(rtracklayer)
library(genomation)

#Setting the right working directory
setwd("/home/jw18713/P1_Full_Run/Datasets/BG3/BG3_histone")

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

# Getting the dm6 genome
genome <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6)

# Choosing a bin size
binsize <- 400

# Creation of the GRange object to contain ChIP metadata
chrs <- seqnames(Dmelanogaster)[1:5]
BG3_tiled <- GRanges()
seqlevels(BG3_tiled) <- chrs
for(i in 1:length(chrs)){
  print(chrs[i])
  start <- seq(1, length(Dmelanogaster[[i]]), by=binsize)
  buffer <- GRanges(seqnames=chrs[i], ranges=IRanges(start = start, end = start + binsize - 1), strand="*")
  BG3_tiled <- c(BG3_tiled, buffer)
}

save(BG3_tiled, file = "/home/jw18713/Archive_Year1/empty_400bp_background.Rda")

background_sample <- sample(BG3_tiled, size = 10000, replace = F)

save(background_sample, file = "/home/jw18713/Archive_Year1/background_sample.Rda")
