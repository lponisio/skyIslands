library(vegan)
library(fields)
library(igraph)
library(betalink)
library(ggplot2)
library(gridExtra)
library(lme4)
library(lmerTest)
library(effects)
library(tidyverse)

source('src/misc.R')
load('../../data/spec.Rdata')


load(file=sprintf("../../data/nets%s%s.Rdata", net.type,
                  paste(species, collapse="")
                  ))


load(file=sprintf('../../data/splev%s%s.Rdata', net.type,
                  paste(species, collapse="")
                  ))

##distance dissimilarity
geo <- unique(spec[, c("Site", "Lat", "Long")])
geo <- geo[!duplicated(geo$Site),]

geo.dist <- rdist.earth(cbind(geo$Long, geo$Lat),
                        cbind(geo$Long, geo$Lat))
colnames(geo.dist) <- rownames(geo.dist) <- geo$Site



plants <- unique(spec$PlantGenusSpecies)
pols <- unique(spec$GenusSpecies)
parasites <- c("AspergillusSpp", "AscosphaeraSpp", "ApicystisSpp",
               "CrithidiaExpoeki", "CrithidiaBombi", "NosemaBombi",
               "NosemaCeranae")

if(species[2] == "Pollinator"){
    lower.level  <- plants
    higher.level <- pols
} else if(species[2] == "Parasite"){
    lower.level  <- pols
    higher.level <- parasites
}
if(net.type == "YrSR"){
    nets.by.SR  <- TRUE
} else {
    nets.by.SR  <- FALSE
}


species.roles <- sp.lev
