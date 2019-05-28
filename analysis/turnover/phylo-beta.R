## setwd('~/Dropbox/skyIslands/')
setwd('analysis/turnover')
source('src/initialize.R')
source('src/phyloIntBeta.R')
library(linkcomm)
library(picante)

## ************************************************************
## prepare link community in terminal
## ************************************************************
## edges.com <- aggregate(list(abund=spec$GenusSpecies),
##                        list(GenusSpecies=spec$GenusSpecies,
##                             PlantGenusSpecies=spec$PlantGenusSpecies),
##                        length)

## lc <- getLinkCommunities(edges.com,
##                          hcmethod = "average",
##                          bipartite=TRUE)
## save(lc, file="saved/lc.Rdata")

## ************************************************************
## turnover of phylo interactions through time
## ************************************************************
load(file="saved/lc.Rdata")
spec$Int <- paste(spec$GenusSpecies,
                  spec$PlantGenusSpecies)

phylo.int.time <- calcCommDis(spec, "Int", lc,
                              abund.w=FALSE,
                              within="Site",
                              between="Year")


phylo.int.space <- calcCommDis(spec, "Int", lc,
                               abund.w=FALSE,
                              within="Year",
                              between="Site",
                              geo.dist=dist.site)
save(phylo.int, file="saved/phyloInt.Rdata")

plot(phylo.int.space$phylo.int$PhyloInt~
         phylo.int.space$phylo.int$Dist)


plot(phylo.int.time$phylo.int$PhyloInt~
         phylo.int.time$phylo.int$Dist)
