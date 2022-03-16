
source("~/P1_Full_Run/Scripts/Model_cleanup/Region_Growth_Function.R")

# Ok, BG3 0.65 first

load("~/P1_Full_Run/Model_outputs/BG3build_BG3glist_0.65.Rda") # glist

no_check <- unlist(glist[lengths(glist) == 1])
to_check <- glist[lengths(glist) > 1]

old <- Sys.time()
grown_regions <- GRanges()
for (i in seq_along(to_check)){
  grown <- Region_Growth(to_check[i], threshold = 0.65, max.gap = 100)
  grown_regions <- c(grown_regions, grown)
}
new <- Sys.time()
print(new - old)

BG3_t65_BG3_trained <- c(no_check, grown_regions)
save(BG3_t65_BG3_trained,
  file = "~/P1_Full_Run/Model_outputs/Thresh_outputs/Threshold_0.65/BG3_t65_BG3_trained.Rda")

rm(list = ls())
source("~/P1_Full_Run/Scripts/Model_cleanup/Region_Growth_Function.R")

# The 0.65 one worked so here's 0.8 for BG3 trained on BG3

load("~/P1_Full_Run/Model_outputs/BG3build_BG3glist_0.8.Rda") # glist

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

BG3_t80_BG3_trained <- c(no_check, grown_regions)
save(BG3_t80_BG3_trained,
  file = "~/P1_Full_Run/Model_outputs/Thresh_outputs/Threshold_0.8/BG3_t80_BG3_trained.Rda")

rm(list = ls())
source("~/P1_Full_Run/Scripts/Model_cleanup/Region_Growth_Function.R")


# From below here I need to do the analysis on S2 trained on BG3_XAI
# 0.65 first


load("~/P1_Full_Run/Model_outputs/BG3build_S2glist_0.65.Rda") # glist

no_check <- unlist(glist[lengths(glist) == 1])
to_check <- glist[lengths(glist) > 1]

old <- Sys.time()
grown_regions <- GRanges()
for (i in seq_along(to_check)){
  grown <- Region_Growth(to_check[i], threshold = 0.65, max.gap = 100)
  grown_regions <- c(grown_regions, grown)
}
new <- Sys.time()
print(new - old)

S2_t65_BG3_trained <- c(no_check, grown_regions)
save(S2_t65_BG3_trained,
  file = "~/P1_Full_Run/Model_outputs/Thresh_outputs/Threshold_0.65/S2_t65_BG3_trained.Rda")

rm(list = ls())
source("~/P1_Full_Run/Scripts/Model_cleanup/Region_Growth_Function.R")


# And S2 for 0.8

load("~/P1_Full_Run/Model_outputs/BG3build_S2glist_0.8.Rda") # glist

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

S2_t80_BG3_trained <- c(no_check, grown_regions)
save(S2_t80_BG3_trained,
  file = "~/P1_Full_Run/Model_outputs/Thresh_outputs/Threshold_0.8/S2_t80_BG3_trained.Rda")

rm(list = ls())
