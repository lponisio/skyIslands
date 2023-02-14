library(vegan)
library(fields)
library(igraph)
library(betalink)
library(ggplot2)
library(gridExtra)
library(bipartite)
## library(effects)

source('src/misc.R')
load('../../data/spec.Rdata')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 0){
    net.type <- args[1]
    species <- args[2]
    sp.level <- args[3]
} else{
    net.type <- "YrSR"
    species <- "Plant"
    sp.level <- "lower.level"
}

plants <- unique(spec$PlantGenusSpecies)
pols <- unique(spec$GenusSpecies)
parasites <- c("AscosphaeraSpp", "ApicystisSpp",
               "CrithidiaExpoeki", "CrithidiaBombi", "NosemaBombi",
               "NosemaCeranae")

if(species == "Plant"){
    species <- c("Plant", "Pollinator")
    lower.level  <- plants
    higher.level <- pols
} else if(species == "Parasite"){
    species <- c("Pollinator", "Parasite")
    lower.level  <- pols
    higher.level <- parasites
}
if(net.type == "YrSR"){
    nets.by.SR  <- TRUE
} else {
    nets.by.SR  <- FALSE
}


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


species.roles <- sp.lev
