# All of the functions I use for constructing tiled genomes
# I need these libraries for the below functions to work
library(GenomicRanges)
library(rtracklayer)
library(genomation)
library(DMRcaller)

# This constructs continuous empty tiled genomes of specific binsizes
# I mainly use this for creating full training datasets
# If single_point == T, constructs single point genomes with a specific gap
# width as defined by binsize
# This single point function is important for DMRcaller to make the ChIP data
# fill happen within my lifetime

empty_tiled_genome <- function(organism, binsize, chr_range, single_point = F){
  chrs <- seqnames(organism)[chr_range]
  tiled_genome <- GRanges()
  seqlevels(tiled_genome) <- chrs
  for(i in 1:length(chrs)){
    print(chrs[i])
    start <- seq(1, length(organism[[i]]), by=binsize)
    if (single_point == F){
      buffer <- GRanges(seqnames=chrs[i],
        ranges=IRanges(start = start, end = start + binsize - 1),
        strand = "*")
    } else if (single_point == T){
      buffer <- GRanges(seqnames=chrs[i],
        ranges=IRanges(start = start, end = start),
        strand = "*")
    }
    tiled_genome <- c(tiled_genome, buffer)
  }
  return(tiled_genome)
}

# Using DMRcaller to average ChIP data over a given region
# When doing ATAC seq, there is score and signalValue
track_calculator <- function(tiled, track, window = 100, target = "score"){
  whole <- c()
  for (c in seqlevels(tiled)){
    chr_tiled <- tiled[seqnames(tiled) == c]
    chr_ovr <- findOverlaps(chr_tiled, track)
    chr_tiled$readsM <- 0
    if (target == "score"){
      chr_tiled$readsM[queryHits(chr_ovr)] <- track$score[subjectHits(chr_ovr)]
    } else if (target == "signalValue"){
      chr_tiled$readsM[queryHits(chr_ovr)] <- track$signalValue[subjectHits(chr_ovr)]
    }
    chr_tiled$readsN <- 1
    chr_tiled$context <- "CG"
    chr_tiled$trinucleotide_context <- "CGG"
    chr_whole <- reduce(chr_tiled, min.gapwidth = window * 2)
    chr_100bp <- computeMethylationProfile(chr_tiled, chr_whole,
      windowSize = window, context = "CG")
    chr_scores <- chr_100bp$Proportion
    chr_scores <- c(chr_scores, 0)
    whole <- c(whole, chr_scores)
  }
  scores <- whole
  return(scores)
}


# Using DMRcaller the way DMRcaller was originally intended
# Slightly modified functions to handle the single point bisulfite data
# without having to create an entire base pair resolution genome

# Function for structuring the original bisulfite tracks the required way
.bisulfite_preprep <- function(tiled, track){
  restructure <- GRanges(seqnames = seqnames(track), ranges = IRanges(start =
  start(track), end = end(track)))
  restructure$readsM <- track$score
  restructure$readsN <- 1
  restructure$context <- "CG"
  restructure$trinucleotide_context <- "CGG"
  return(restructure)
}

# Bisulfite track calculator function in full
Bisulfite_calculator <- function(tiled, track, win = 100){
  processed_track <- .bisulfite_preprep(tiled, track)
  whole <- c()
  for (c in seqlevels(tiled)){
    chr <- processed_track[seqnames(processed_track) == c]
    chr_whole <- reduce(tiled[seqnames(tiled) == c], min.gapwidth = win * 2)
    chr_100bp <- computeMethylationProfile(chr, chr_whole,
      windowSize = win, context = "CG")
    chr_100bp$Proportion[chr_100bp$Proportion == "NaN"] <- 0
    chr_scores <- chr_100bp$Proportion
    chr_scores <- c(chr_scores, 0)
    whole <- c(whole, chr_scores)
  }
  scores <- whole
  return(scores)
}

# Function for selecting a specific number of regions with non-zero values
f_select <- function(x, n){
  check <- mcols(x)[grep("log2", names(mcols(x)))-1]
  atac_less <- check[,1:length(check)-1]
  greater <- apply(as.matrix(atac_less)>0, 1, sum)
  to_select <- greater >= n
  selected <- x[to_select]

  # Coverage Report
  sel_perc <- round(100/length(x)*length(selected), 2)
  print(paste0("The selected threshold covers ", sel_perc,
    "% of the genome."))

  # STARR report
  starr_perc <- round(100/sum(x$STARR_seq_binary)*sum(selected$STARR_seq_binary), 2)
  print(paste0("The selected threshold covers ", starr_perc,
    "% of positive STARR-seq bins."))

  return(selected)
}
