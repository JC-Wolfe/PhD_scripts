# Get the functions
source("~/Human_Datasets/Source_scripts/Human_glist_functions.R")
source("~/Human_Datasets/Source_scripts/expression_analysis_functon.R")
# Load the H9 objects
# Proximal and Distal is 1
# Distal Only is 2
# Proximal Only is 3
# Neither is 4
# So distal_only will be list[[2]]
# proximal_all will be sort(c(list[[1]],list[[3]]))

setwd("~/Pre_viva/Hi-C")

# Putative
load("contact_Rda/H9_background.Rda")

# loading distal_contacts
load("contact_Rda/distal_contacts.Rda")

# expression data
load("expression_data.Rda")
# And finally, running the things themselves
# Putative Enhancers

# For H9 background enhancers
# Distal
H9_background_dist <- distal_expression(H9_background[[2]], g_exp = expression_data,
  distal_contacts = distal_contacts)

H9_background_distal_EA <- expression_analysis(H9_background_dist)

save(H9_background_distal_EA,
  file = "~/Pre_viva/Hi-C/grange_lists/H9_background_distal_EA.Rda")
# Proximal
H9_background_prox <- proximal_expression(sort(c(H9_background[[1]],H9_background[[3]])),
  g_exp = expression_data)

H9_background_proximal_EA <- expression_analysis(H9_background_prox)

save(H9_background_proximal_EA,
  file = "~/Pre_viva/Hi-C/grange_lists/H9_background_proximal_EA.Rda")
