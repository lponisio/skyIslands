library(igraph)
library(bipartite)
library(lme4)
library(lmerTest)
library(RColorBrewer)

source('src/CalcMetrics.R')
source('src/misc.R')

save.path <- 'saved'

load('../../data/spec.Rdata')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 0){
    net.type <- args[1]
    species <- args[2]
} else{
    net.type <- "YrSR"
    species <- "Parasite"
}


plants <- unique(spec$PlantGenusSpecies)
pols <- unique(spec$GenusSpecies)
parasites <- c("AspergillusSpp", "AscosphaeraSpp", "ApicystisSpp",
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


