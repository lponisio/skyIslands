## This script calculates the pollinator role/network niche
## differences across sites/years
## setwd('~/Dropbox/skyIslands')
rm(list=ls())
setwd('analysis/variability')
source('src/initialize.R')

species.roles <- calcSpec(nets, spec, dist.metric="chao")

## vector of pca loadings of interest
loadings <- c(1)
metrics <- c("rare.degree",
             "weighted.betweenness",
             "weighted.closeness",
             "niche.overlap",
             "species.strength",
             "d")

pol <- species.roles[species.roles$speciesType == "pollinator",]
pol.pca.scores <- calcNetworkPca (species.roles=pol,
                                  metrics= metrics,
                                  loadings=loadings)


pcas.species <- split(pol.pca.scores$pcas,
                      pol.pca.scores$pcas$GenusSpecies)

pcas.species <- lapply(pcas.species, calcDiffPcas, geo.dist)
pcas.species <- pcas.species[!sapply(pcas.species, is.null)]
pcas.species <- do.call(rbind, pcas.species)
rownames(pcas.species) <- NULL


## turnover through space
pcas.diffs.same.year <- pcas.species[apply(pcas.species, 1,
                                           function(x) x["Year1"] ==
                                                       x["Year2"] &
                                                       x["SR1"] ==
                                                       x["SR2"]),]


p <- ggplot(pcas.diffs.same.year,
       aes(x=GeoDist, y=diffPca,
           pch=Year1, col=GenusSpecies)) + geom_point() + geom_line() +
    labs(x="Geographic Distance", y="Network role change") +
    lims(y=c(0,1))

p + theme(legend.position="top", legend.box = "horizontal")

