
enhancers <- read.csv("~/Review_Paper_Plots/arch_prot_rules/data/ChIP_enrichment_95perc_quartile_table_enhancer_width_window_raw.csv")
background <- read.csv("~/Review_Paper_Plots/arch_prot_rules/data/ChIP_enrichment_95perc_quartile_table_2Kb_genome_raw.csv")

enhancers <- enhancers[,c(2:21,30:32)]
background <- background[,c(1:20,29:30)]

background_norm <- function(x, y){
  for (i in seq(1,length(x)){
    mx <- max(y[,i])
    mn <- min(y[,i])
    x[,i] <- (x[,i] - mn)/(mx - mn)
  }
  return(x)
}

norm_enhancers <- enhancers
norm_background <- background

norm_enhancers[2:22] <- background_norm(enhancers[2:22], background[2:22])
norm_background[2:22] <- background_norm(background[2:22], background[2:22])
