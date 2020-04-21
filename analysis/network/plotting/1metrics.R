## setwd('~/Dropbox/skyIslands/')
rm(list=ls())
setwd('analysis/network')
source('src/initialize.R')
source("src/misc.R")
source("plotting/src/predictIntervals.R")
source("plotting/src/CIplotting.R")
source("plotting/src/diagnostics.R")
source("plotting/src/plotNetworkMets.R")

species <- c("Plant", "Pollinator")
## species <- c("Pollinator", "Parasite")

xvar <- c("Lat")
xlabel <- "Latitude"

load(file=file.path(save.path,
                    sprintf('mods/metrics_%s_%s.Rdata',
                            paste(species, collapse=""),
                            xvar)))

## ************************************************************
## network metrics
## ************************************************************

ys <- names(mods.div)

ylabs <- c("Plant \n niche overlap",
           "Pollinatator \n niche overlap",
           "Plant clustering",
           "Pollinator clustering",
           "Plant richness",
           "Pollinator richness",
           "Nestedness (z)",
           "Reciprocal specialization (H2)")
names(ylabs) <- ys


x <- xvar


pdf.f(plotNetworkMets,
      file=file.path("figures",
                     sprintf("%s_%s.pdf", xvar,
                             paste(species, collapse=""))),
      width=8.5, height=11)

