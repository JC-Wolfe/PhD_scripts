
library(GenomicRanges)
library(rtracklayer)

.expression_extraction <- function(gr, expression_data){
  ovr <- findOverlaps(gr, expression_data)
  gr <- gr[queryHits(ovr)]
  gr$contact_chr <- expression_data[subjectHits(ovr)]$contact_chr
  gr$contact_start <- expression_data[subjectHits(ovr)]$contact_start
  gr$contact_end <- expression_data[subjectHits(ovr)]$contact_end
  gr$Gene_short_name <- expression_data[subjectHits(ovr)]$Gene_short_name
  gr$FPKM <- expression_data[subjectHits(ovr)]$FPKM
  gr <- gr[!duplicated(gr$FPKM)]
  return(gr)
}

distal_expression <- function(x, g_exp, distal_contacts){

  exp_proms <- g_exp
  start(exp_proms[strand(exp_proms) == "+"]) <- start(exp_proms[strand(exp_proms) == "+"])-250
  end(exp_proms[strand(exp_proms) == "+"]) <- start(exp_proms[strand(exp_proms) == "+"])
  end(exp_proms[strand(exp_proms) == "-"]) <- end(exp_proms[strand(exp_proms) == "-"])+250
  start(exp_proms[strand(exp_proms) == "-"]) <- end(exp_proms[strand(exp_proms) == "-"])

  # So lets see where we have promoter contacts then
  exp_test <- GRanges(seqnames = distal_contacts$contact_chr,
    ranges = IRanges(distal_contacts$contact_start, distal_contacts$contact_end))
  seqlevelsStyle(exp_test) <- "UCSC"

  exp_over <- findOverlaps(exp_test, exp_proms)
  contacted_expressed <- distal_contacts[queryHits(exp_over)]
  contacted_expressed$Gene_short_name <- exp_proms$Gene_short_name[subjectHits(exp_over)]
  contacted_expressed$FPKM <- exp_proms$FPKM[subjectHits(exp_over)]

  # This is a GRange list that will contain all of the target enhancers position by position
  list_distal_only <- split(x, as.factor(x))

  old <- Sys.time()
  glist <- mclapply(list_distal_only, FUN = .expression_extraction,
    mc.cores = 30, expression_data = contacted_expressed)
  new <- Sys.time()
  print(old - new)
  return(glist)
}

# Proximal functions are down here

proximal_expression <- function(x, g_exp){

  # So lets see where we have promoter contacts then
  enh_reach <- GRanges(seqnames = seqnames(x), ranges = IRanges(start(x)-5000, end(x)+5000))
  exp_proms <- g_exp
  start(exp_proms[strand(exp_proms) == "+"]) <- start(exp_proms[strand(exp_proms) == "+"])-250
  end(exp_proms[strand(exp_proms) == "+"]) <- start(exp_proms[strand(exp_proms) == "+"])
  end(exp_proms[strand(exp_proms) == "-"]) <- end(exp_proms[strand(exp_proms) == "-"])+250
  start(exp_proms[strand(exp_proms) == "-"]) <- end(exp_proms[strand(exp_proms) == "-"])
  exp_over <- findOverlaps(enh_reach, exp_proms)
  contacted_expressed <- enh_reach[queryHits(exp_over)]
  contacted_expressed$Gene_short_name <- exp_proms$Gene_short_name[subjectHits(exp_over)]
  contacted_expressed$FPKM <- exp_proms$FPKM[subjectHits(exp_over)]

  # This is a GRange list that will contain all of the target enhancers position by position
  list_proximal <- split(x, as.factor(x))

  old <- Sys.time()
  glist <- mclapply(list_proximal, FUN = .expression_extraction,
    mc.cores = 30, expression_data = contacted_expressed)
  new <- Sys.time()
  print(old - new)
  return(glist)
}
