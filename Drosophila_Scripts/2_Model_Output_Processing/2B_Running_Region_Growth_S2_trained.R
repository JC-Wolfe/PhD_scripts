
source("~/P1_Full_Run/Scripts/Model_cleanup/Region_Growth_Function.R")

# Ok, BG3 0.75 first

load("~/P1_Full_Run/S2_Model_outputs/S2build_BG3glist_0.75.Rda") # glist

no_check <- unlist(glist[lengths(glist) == 1])
to_check <- glist[lengths(glist) > 1]

old <- Sys.time()
grown_regions <- GRanges()
for (i in seq_along(to_check)){
  grown <- Region_Growth(to_check[i], threshold = 0.75, max.gap = 100)
  grown_regions <- c(grown_regions, grown)
}
new <- Sys.time()
print(new - old)

BG3_t75_S2_trained <- c(no_check, grown_regions)
save(BG3_t75_S2_trained,
  file = "~/P1_Full_Run/S2_Model_outputs/Thresh_outputs/Threshold_0.75/BG3_t75_S2_trained.Rda")

rm(list = ls())
source("~/P1_Full_Run/Scripts/Model_cleanup/Region_Growth_Function.R")

# The 0.75 one worked so here's 0.8 for BG3 trained on BG3

load("~/P1_Full_Run/S2_Model_outputs/S2build_BG3glist_0.8.Rda") # glist

no_check <- unlist(glist[lengths(glist) == 1])
to_check <- glist[lengths(glist) > 1]

old <- Sys.time()
grown_regions <- GRanges()
for (i in seq_along(to_check)){
  grown <- Region_Growth(to_check[i], threshold = 0.8, max.gap = 100)
  grown_regions <- c(grown_regions, grown)
}
new <- Sys.time()
print(new - old)

BG3_t80_S2_trained <- c(no_check, grown_regions)
save(BG3_t80_S2_trained,
  file = "~/P1_Full_Run/S2_Model_outputs/Thresh_outputs/Threshold_0.8/BG3_t80_S2_trained.Rda")

rm(list = ls())
source("~/P1_Full_Run/Scripts/Model_cleanup/Region_Growth_Function.R")


# From below here I need to do the analysis on S2 trained on BG3_XAI
# 0.75 first


load("~/P1_Full_Run/S2_Model_outputs/S2build_S2glist_0.75.Rda") # glist

no_check <- unlist(glist[lengths(glist) == 1])
to_check <- glist[lengths(glist) > 1]

old <- Sys.time()
grown_regions <- GRanges()
for (i in seq_along(to_check)){
  grown <- Region_Growth(to_check[i], threshold = 0.75, max.gap = 100)
  grown_regions <- c(grown_regions, grown)
}
new <- Sys.time()
print(new - old)

S2_t75_S2_trained <- c(no_check, grown_regions)
save(S2_t75_S2_trained,
  file = "~/P1_Full_Run/S2_Model_outputs/Thresh_outputs/Threshold_0.75/S2_t75_S2_trained.Rda")

rm(list = ls())
source("~/P1_Full_Run/Scripts/Model_cleanup/Region_Growth_Function.R")


# And S2 for 0.8

load("~/P1_Full_Run/S2_Model_outputs/S2build_S2glist_0.8.Rda") # glist

no_check <- unlist(glist[lengths(glist) == 1])
to_check <- glist[lengths(glist) > 1]

old <- Sys.time()
grown_regions <- GRanges()
for (i in seq_along(to_check)){
  grown <- Region_Growth(to_check[i], threshold = 0.8, max.gap = 100)
  grown_regions <- c(grown_regions, grown)
}
new <- Sys.time()
print(new - old)

S2_t80_S2_trained <- c(no_check, grown_regions)
save(S2_t80_S2_trained,
  file = "~/P1_Full_Run/S2_Model_outputs/Thresh_outputs/Threshold_0.8/S2_t80_S2_trained.Rda")

rm(list = ls())
