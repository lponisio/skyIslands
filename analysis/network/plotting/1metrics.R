## setwd('~/Dropbox/skyIslands/')
rm(list=ls())
setwd('analysis/network')
source('src/initialize.R')
source("src/misc.R")
source("plotting/src/predictIntervals.R")
source("plotting/src/CIplotting.R")
source("plotting/src/diagnostics.R")
source("plotting/src/plotNetworkMets.R")

## species <- c("Plant", "Pollinator")
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


x <- xvar

if(species[1] == "Plant"){

    ylabs <- c("Plant \n niche overlap",
               "Pollinatator \n niche overlap",
               "Plant clustering",
               "Pollinator clustering",
               "Mean partners plants",
               "Mean parnters pollinators",
               "Plant richness",
               "Pollinator richness",
               "Nestedness (zNODF)",
               "Reciprocal specialization (zH2')")
    names(ylabs) <- ys

    yrs <- c("2012", "2017", "2018")
    years <- list(niche.overlap.LL=yrs[1],
                  niche.overlap.HL=yrs[1],
                  weighted.cluster.coefficient.LL=yrs,
                  weighted.cluster.coefficient.HL=yrs,
                  mean.number.of.links.LL=yrs,
                  mean.number.of.links.HL=yrs,
                  number.of.species.LL=yrs,
                  number.of.species.HL=yrs,
                  zweighted.NODF=yrs,
                  H2=yrs)

    pdf.f(plotNetworkMetsPlantPol,
          file=file.path("figures",
                         sprintf("%s_%s.pdf", xvar,
                                 paste(species, collapse=""))),
          width=8.5, height=11)


    pdf.f(plotNetworkMetsPlantPolDOY,
          file=file.path("figures",
                         sprintf("%s_%s.pdf", "Doy",
                                 paste(species, collapse=""))),
          width=8.5, height=11)

} else if(species[2] == "Parasite"){
    ylabs <- c("Pollinator \n niche overlap",
               "Parasite \n niche overlap",
               "Pollinator clustering",
               "Parasite clustering",
               "Mean partners pollinators",
               "Mean parnters parasites",
               "Pollinator richness",
               "Parasite richness",
               "Nestedness (zNODF)",
               "Reciprocal specialization (zH2')")
    names(ylabs) <- ys

    yrs <- c("2018")
    years <- list(niche.overlap.LL=yrs,
                  niche.overlap.HL=yrs,
                  weighted.cluster.coefficient.LL=yrs,
                  weighted.cluster.coefficient.HL=yrs,
                  mean.number.of.links.LL=yrs,
                  mean.number.of.links.HL=yrs,
                  number.of.species.LL=yrs,
                  number.of.species.HL=yrs,
                  zweighted.NODF=yrs,
                  H2=yrs)

    pdf.f(plotNetworkMetsPlantPol,
          file=file.path("figures",
                         sprintf("%s_%s.pdf", xvar,
                                 paste(species, collapse=""))),
          width=8.5, height=11)


    pdf.f(plotNetworkMetsPlantPolDOY,
          file=file.path("figures",
                         sprintf("%s_%s.pdf", "Doy",
                                 paste(species, collapse=""))),
          width=8.5, height=11)

}
