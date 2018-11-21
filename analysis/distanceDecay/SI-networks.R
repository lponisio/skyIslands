rm(list=ls())
library(vegan)
library(fields)
library(igraph)
library(linkcomm)
setwd("~/Dropbox/skyIslands/analysis/distanceDecay")
source('src/misc.R')
source('src/plot-mod.R')
source('~/Dropbox/network_assembly/simulation/src/all/CalcMetrics.R')

spec <-
  read.csv("../data/spec.csv")

prep.comm <- aggregate(spec$GenSp,
                       list(site=spec$PlantGenSp,
                            sp=spec$GenSp), length)

comm <- samp2site.spp(prep.comm$site, prep.comm$sp, prep.comm$x)

g <- mut.adj(comm)
weights <- as.vector(comm)
weights <- weights[weights != 0]
mod <- edge.betweenness.community(g, weights=weights, directed=FALSE)

## ************************************************************
## prepare link community in terminal
## ************************************************************

edges.com <-cbind(as.character(spec$GenSp),
                  as.character(spec$PlantGenSp))

lc <- getLinkCommunities(edges.com, hcmethod = "average",
                         bipartite=TRUE)
save(lc, file="src/lc.Rdata")

## ************************************************************

load("src/lc.Rdata")

nested.comm <- getAllNestedComm(lc)
cc <- getCommunityCentrality(lc)
cm <- getCommunityConnectedness(lc, conn = "mod")

path <- '../figures/modularity'

plot.all()
