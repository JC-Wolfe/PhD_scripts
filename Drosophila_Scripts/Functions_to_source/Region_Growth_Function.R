library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)

.sumReadsConfPerc1 <- function(methylationData){   #Function to compute sums within below function
  return(sum(methylationData$sumConfPerc1))
}

.ConfPercAvg <- function(methylationData, regions){

  overlaps <- findOverlaps(methylationData, regions, ignore.strand = TRUE)
  methylationDataContextList <- IRanges::splitAsList(methylationData[queryHits(overlaps)],  subjectHits(overlaps))
  regionsIndexes <- as.integer(names(methylationDataContextList))

  regions$sumConfPerc1 <- rep(0, times=length(regions)) #Confidence percentage 1

  if(length(regionsIndexes) > 0){
    regions$sumConfPerc1[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsConfPerc1) #Confidence percentage 1
  }

  regions$avgConfPerc1 <- regions$sumConfPerc1 / (width(regions)/10)

  return(regions)
}

Region_Growth <- function(x, threshold, max.gap = 100){

  x <- unlist(x) # This is to get the object out of the glist (I don't need this if I can't sapply)
  n <- length(x) # This is the length number for the recursive function

  # First I have to check if the entire region can be combined within the threshold
  all_reduced <- .ConfPercAvg(x, reduce(x, min.gapwidth = max.gap + 1))

  # Now I can stop here for speed if all regions can be combined
  if (all_reduced$avgConfPerc1 >= threshold){
    return(all_reduced)
  }

  # This is what happens if not every region can be combined
  else{

    all_results <- x # A place to store my results (contains uncombined regions)
    n1 <- n # Starting from n-1 since not all can be combined

    # Gradually decreasing slice size until a minimum of 2 is reached
    while(n > 2){
      n <- n-1 # Reducing the length of the slices
      start_vector <- seq(1,(n1-n+1),by=1) # Getting the start of each slice
      end_vector <- seq(n1-(n1-n),n1,by=1) # Getting the end of each slice
      reduced_candidates <- GRanges() # Holding place for reduced region stats

      # Slices, for example if we're looking at 5 we'd have 1-4 and 2-5 at n = 4,
      # then, as n decreases we'd have 1-3, 2-4, and 3-5 at n = 3, etc.
      for (i in start_vector){
        candidate <- x[i: end_vector[i]] # Individual slices
        reduced <- .ConfPercAvg(candidate, reduce(candidate, min.gapwidth = max.gap + 1))
        reduced_candidates <- c(reduced_candidates, reduced)
      }
      all_results <- c(all_results, reduced_candidates) # Appending all results
    }
    above_threshold <- all_results[all_results$avgConfPerc1 >= threshold]
    reduced_regions <- .ConfPercAvg(above_threshold, reduce(above_threshold))
    return(reduced_regions)
  }
}
