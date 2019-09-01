## setwd('~/Dropbox/skyislands')
rm(list=ls())
setwd('analysis/variability')
source('src/initialize.R')
type <- "all"

species.roles <- calcSpec(nets, spec, dist.metric="chao")

## vector of pca loadings of interest
loadings <- c(1)
metrics <- c("rare.degree", "weighted.betweenness", "weighted.closeness",
             "niche.overlap", "d")

## the metrics used in the PCA
var.method <- cv
ave.method <- mean

## *********************************************************************
## pollinators
## *********************************************************************

pol <- species.roles[species.roles$speciesType == "pollinator",]
pol.pca.scores <- calcPcaMeanVar(species.roles=pol,
                                 var.method=var.method,
                                 ave.method=ave.method,
                                 metrics= metrics,
                                 loadings=loadings,
                                 agg.col = "Year")

autoplot(pol.pca.scores$pca.loadings$'2017', loadings=TRUE,
         ## loadings.label = TRUE,
         ## loadings.label.size = 2,
         loadings.colour = 'blue')


## *********************************************************************
## plants
## *********************************************************************

plants <- species.roles[species.roles$speciesType == "plant",]
plant.pca.scores <- calcPcaMeanVar(species.roles=plants,
                                   var.method=var.method,
                                   metrics= metrics,
                                   loadings=loadings,
                                   ave.method=ave.method,
                                   agg.col = "Year")

save(pol.pca.scores, plant.pca.scores, file="saved/results/pcaVar.Rdata")
