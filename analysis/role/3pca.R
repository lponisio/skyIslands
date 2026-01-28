## This script calculates the pollinator role/network niche
## variability between sites within a year.
rm(list=ls())
library(ggfortify)
library(bipartite)
library(fossil)
setwd("C:/")
source("lab_paths.R")
local.path

dir.bombus <- file.path(local.path, "skyIslands")
setwd(dir.bombus)
getwd()

this.script <- "role"
source('analysis/role/src/initialize.R')
type <- "all"

## drop empty cells and NA
spec.net <-  spec.net[spec.net$GenusSpecies != '',]
spec.net <-  spec.net[spec.net$PlantSpecies != '',]
spec.net <-  spec.net[spec.net$PlantGenus != '',]
spec.net <- spec.net[!is.na(spec.net$GenusSpecies),]
spec.net <- spec.net[!is.na(spec.net$PlantSpecies),]
spec.net <- spec.net[!is.na(spec.net$PlantGenus),]

any(spec.net$GenusSpecies == "")
any(is.na(spec.net$GenusSpecies))
any(spec.net$PlantGenus == "")
any(is.na(spec.net$PlantSpecies))


nets <- nets[ names(nets) != "PL.2017" ]

spec.net$Year <- as.numeric(spec.net$Year)
sum(is.na(spec.net$Year))  # should be 0


names(nets)


## calculate species roles
species.roles <- calcSpec(nets, spec.net, dist.metric="chao")
species.roles.tidy <- species.roles[species.roles$speciesType == "plant",]


str(spec.net)

## vector of pca loadings of interest
loadings <- c(1)
metrics <- c("rare.degree",
             "weighted.betweenness",
             "weighted.closeness",
             "niche.overlap",
             "species.strength",
             "d")

## Metrics used in the PCA
var.method <- cv
ave.method <- mean

## PCA 
plant.pca.scores <- calcPcaMeanVar(species.roles=plant, 
                                 var.method=var.method,
                                 ave.method=ave.method,
                                 metrics= metrics,
                                 loadings=loadings,
                                 agg.col = "Year")


autoplot(plant.pca.scores$'2019'$pca.loadings, loadings=TRUE,
         loadings.colour = 'blue')


save(plant.pca.scores,  file="saved/results/pcaVar.Rdata")
