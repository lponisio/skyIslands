## setwd("~/Dropbox/skyIslands/")
rm(list=ls())
setwd("analysis/turnover")
source("src/initialize.R")
source("src/chao.R")
source("src/betalink.R")


nets.graph <- prepare_networks(nets)


beta.net <- network_betadiversity(nets.graph)

beta.net$propOS <- beta.net$OS/beta.net$WN
beta.net$propST <- beta.net$ST/beta.net$WN
