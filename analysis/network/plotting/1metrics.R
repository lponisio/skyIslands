## setwd('~/Dropbox/skyIslands/')
rm(list=ls())
setwd('analysis/network')
source('src/initialize.R')
source("src/misc.R")
source("plotting/src/predictIntervals.R")
source("plotting/src/CIplotting.R")
source("plotting/src/diagnostics.R")
source("plotting/src/plotNetworkMets.R")

xvars <- c("Lat")
xlabel <- "Latitude"

## ************************************************************
## network metrics
## ************************************************************

load('saved/mods/metrics.Rdata')

ys <- names(mods.div)
ylabs <- c("Plant \n complementarity",
           "Pollinatator \n complementarity",
           "Reciprocal \n specialization",
           "Niche breadth",
           "Plant \n niche overlap", "Pollinator \n niche overlap",
           "Plant \n richness", "Pollinator \n richness")
names(ylabs) <- ys

x <- xvars


pdf.f(plotNetworkMets, file=file.path("figures",
                                      sprintf("%s.pdf", xlabel)),
      width=8.5, height=11)

