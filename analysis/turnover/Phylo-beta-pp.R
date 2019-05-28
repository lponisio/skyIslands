## setwd("~/Dropbox/skyIslands")
setwd("analysis/turnover")
rm(list=ls())
library(vegan)
library(fields)
library(linkcomm)
library(picante)
source('src/misc.R')
source('src/Phylo-beta-int.R')
source('src/Phylo-beta-sp.R')
source('src/misc.R')

geo <-
  read.csv("../../data/relational/data/relational/tables/geography.csv")

spec <-
  read.csv("../data/spec.csv")

## ************************************************************
## prepare link community in terminal
## ************************************************************

edges.com <-cbind(as.character(spec$GenusSpecies),
                  as.character(spec$PlantGenusSpecies))

lc <- getLinkCommunities(edges.com, hcmethod = "average",
                         bipartite=TRUE)
save(lc, file="src/lc.Rdata")

## ************************************************************

## distance dissimilarity
dist.site <- rdist.earth(cbind(geo$Long, geo$Lat),
                         cbind(geo$Long,
                               geo$Lat), miles=FALSE)

c.dist <- dist.site[lower.tri(dist.site)]

## dissimilarity of plants, pol, int

int <- int.dis(spec, c.dist, abund.w=TRUE,
               path ='../figures/distanceDecay')

sp <- pp.dis(spec, types= c("GenSp", "PlantGenSp", "Int"),
             c.dist=c.dist, sub="all",
             path ='../figures/distanceDecay')

