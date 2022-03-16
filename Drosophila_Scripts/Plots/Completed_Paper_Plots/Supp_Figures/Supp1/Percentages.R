load("BG3_putative_list.Rda")
load("BG3_shared_list.Rda")
load("BG3_starr_res_list.Rda")
load("BG3_bg_list.Rda")
load("S2_putative_list.Rda")
load("S2_shared_list.Rda")
load("S2_starr_res_list.Rda")
load("S2_bg_list.Rda")

#percentage getting

BG3_putative_total <- length(unlist(BG3_putative))
BG3_putative_distal <- length(BG3_putative[[2]])
BG3_putative_prox <- length(c(BG3_putative[[1]], BG3_putative[[3]]))
BG3_putative_neither <- length(BG3_putative[[4]])
BG3_putative_prox_percent <- 100/BG3_putative_total*BG3_putative_prox
BG3_putative_dist_percent <- 100/BG3_putative_total*BG3_putative_distal
BG3_putative_neither_percent <- 100/BG3_putative_total*BG3_putative_neither


BG3_common_total <- length(unlist(BG3_shared))
BG3_common_distal <- length(BG3_shared[[2]])
BG3_common_prox <- length(c(BG3_shared[[1]], BG3_shared[[3]]))
BG3_common_neither <- length(BG3_shared[[4]])
BG3_common_prox_percent <- 100/BG3_common_total*BG3_common_prox
BG3_common_dist_percent <- 100/BG3_common_total*BG3_common_distal
BG3_common_neither_percent <- 100/BG3_common_total*BG3_common_neither

S2_putative_total <- length(unlist(S2_putative))
S2_putative_distal <- length(S2_putative[[2]])
S2_putative_prox <- length(c(S2_putative[[1]], S2_putative[[3]]))
S2_putative_neither <- length(S2_putative[[4]])
S2_putative_prox_percent <- 100/S2_putative_total*S2_putative_prox
S2_putative_dist_percent <- 100/S2_putative_total*S2_putative_distal
S2_putative_neither_percent <- 100/S2_putative_total*S2_putative_neither

S2_common_total <- length(unlist(S2_shared))
S2_common_distal <- length(S2_shared[[2]])
S2_common_prox <- length(c(S2_shared[[1]], S2_shared[[3]]))
S2_common_neither <- length(S2_shared[[4]])
S2_common_prox_percent <- 100/S2_common_total*S2_common_prox
S2_common_dist_percent <- 100/S2_common_total*S2_common_distal
S2_common_neither_percent <- 100/S2_common_total*S2_common_neither

# Background stuff

BG3_bg_total <- length(unlist(BG3_bg))
BG3_bg_distal <- length(BG3_bg[[2]])
BG3_bg_prox <- length(c(BG3_bg[[1]], BG3_bg[[3]]))
BG3_bg_neither <- length(BG3_bg[[4]])
BG3_bg_prox_percent <- 100/BG3_bg_total*BG3_bg_prox
BG3_bg_dist_percent <- 100/BG3_bg_total*BG3_bg_distal
BG3_bg_neither_percent <- 100/BG3_bg_total*BG3_bg_neither

S2_bg_total <- length(unlist(S2_bg))
S2_bg_distal <- length(S2_bg[[2]])
S2_bg_prox <- length(c(S2_bg[[1]], S2_bg[[3]]))
S2_bg_neither <- length(S2_bg[[4]])
S2_bg_prox_percent <- 100/S2_bg_total*S2_bg_prox
S2_bg_dist_percent <- 100/S2_bg_total*S2_bg_distal
S2_bg_neither_percent <- 100/S2_bg_total*S2_bg_neither

# STARR-seq

BG3_starr_res_total <- length(unlist(BG3_starr_res))
BG3_starr_res_distal <- length(BG3_starr_res[[2]])
BG3_starr_res_prox <- length(c(BG3_starr_res[[1]], BG3_starr_res[[3]]))
BG3_starr_res_neither <- length(BG3_starr_res[[4]])
BG3_starr_res_prox_percent <- 100/BG3_starr_res_total*BG3_starr_res_prox
BG3_starr_res_dist_percent <- 100/BG3_starr_res_total*BG3_starr_res_distal
BG3_starr_res_neither_percent <- 100/BG3_starr_res_total*BG3_starr_res_neither

S2_starr_res_total <- length(unlist(S2_starr_res))
S2_starr_res_distal <- length(S2_starr_res[[2]])
S2_starr_res_prox <- length(c(S2_starr_res[[1]], S2_starr_res[[3]]))
S2_starr_res_neither <- length(S2_starr_res[[4]])
S2_starr_res_prox_percent <- 100/S2_starr_res_total*S2_starr_res_prox
S2_starr_res_dist_percent <- 100/S2_starr_res_total*S2_starr_res_distal
S2_starr_res_neither_percent <- 100/S2_starr_res_total*S2_starr_res_neither
