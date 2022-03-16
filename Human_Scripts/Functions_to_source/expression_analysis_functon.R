
.maxExp <- function(x){
  return(max(x$FPKM, na.rm = T))
}

.meanExp <- function(x){
  return(mean(x$FPKM, na.rm = T))
}

.medExp <- function(x){
  return(median(x$FPKM, na.rm = T))
}

expression_analysis <- function(x, cores = 30){
  # Creating an object to fill with expression information for later plots
  e_ob <- list()

  # Finding which regions have at least one promoter contact
  contacts <- x[!lengths(x) == 0]
  no_contacts <- x[lengths(x) == 0]

  # Adding this information to the list to be output
  e_ob$count_with_contacts <- length(contacts)
  e_ob$count_without_contacts <- length(no_contacts)
  e_ob$contact_proportion <- length(contacts)/(length(contacts) + length(no_contacts))

  # Finding the max, mean and median expression vectors
  e_ob$max_Expression <- unlist(mclapply(contacts, FUN = .maxExp,
    mc.cores = cores))
  e_ob$mean_Expression <- unlist(mclapply(contacts, FUN = .meanExp,
    mc.cores = cores))
  e_ob$median_Expression <- unlist(mclapply(contacts, FUN = .medExp,
    mc.cores = cores))

  return(e_ob)
}
