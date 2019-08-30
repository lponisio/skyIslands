setwd('~/Dropbox/skyIslands/')
rm(list=ls())
setwd('analysis/network')
source('src/initialize.R')
source("src/misc.R")
source("plotting/src/predictIntervals.R")
source("plotting/src/CIplotting.R")
source("plotting/src/diagnostics.R")
source("plotting/src/plotNetworkMets.R")

xvars <- c("Lat")
type <- "all"
xlabel <- "Latitude"

## ************************************************************
## network metrics
## ************************************************************

load('saved/mods/metrics.Rdata')

ys <- names(mods.div)
ylabs <- ys
names(ylabs) <- ys
x <- xvars
## x2 <- xvars[2]

pdf.f(plotNetworkMets, file=file.path("figures",
                                      sprintf("%s.pdf", type)),
      width=5, height=8)
