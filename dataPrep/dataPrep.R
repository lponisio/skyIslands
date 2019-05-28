rm(list=ls())
library(igraph)
library(vegan)
library(fields)
library(bipartite)
setwd("~/Dropbox/skyIslands/dataPrep")
source('relational/relational_prep.R')
setwd("~/Dropbox/skyIslands/dataPrep")
source('relational/relational_make.R')
setwd("~/Dropbox/skyIslands/dataPrep")
source('relational/make_traditional.R')

setwd("~/Dropbox/skyIslands/dataPrep")
src.dir <- '../../skyIslands_saved/data/relational/relational/traditional/'
spec <-
    read.csv(file.path(src.dir, "specimens-complete.csv"))

source("src/misc.R")
source("src/prepNets.R")
source("src/specialization.R")

## did not complete full sampling rounds in any of these sites. Was
## just scouting.
## can keep UK and SS when more species are IDed
site.2.drop <- c("JM", "CC", "UK", "SS")
spec <- spec[!spec$Site %in% site.2.drop,]
spec <- droplevels(spec)

## get specimen data ready
spec$GenusSpecies <- fix.white.space(paste(spec$Genus,
                          spec$Species,
                          spec$SubSpecies))

spec$PlantGenusSpecies <-  fix.white.space(paste(spec$PlantGenus,
                                          spec$PlantSpecies,
                                          spec$PlantVar,
                                          spec$PlantSubSpecies))

spec$Int <-  fix.white.space(paste(spec$GenusSpecies,
                                   spec$PlantGenusSpecies))
spec$IntGen <-  fix.white.space(paste(spec$Genus,
                                      spec$PlantGenus))

spec$Date <- as.Date(spec$Date, format='%m/%d/%y')
spec$Doy <- as.numeric(strftime(spec$Date, format='%j'))
spec$Year <- as.numeric(format(spec$Date,'%Y'))

## drop non-bees
spec <- spec[spec$Family %in% c("Andrenidae", "Apidae",
                                "Colletidae", "Halictidae",
                                "Megachilidae"),]


## for networks, drop specimens withoutplant IDs
spec <- spec[spec$PlantGenusSpecies != "",]
spec <- spec[spec$GenusSpecies != "",]

## drop the the 2017 sample of PL because it was on fire for other
## sampling rounds and there was basically nothing blooming the first
## round
spec <- spec[!(spec$Site == "PL" & spec$Year == "2017"),]

save(spec, file="../data/spec.Rdata")


### networks
spec$YearSR <- paste(spec$Year, spec$SampleRound, sep=".")

nets <- breakNet(spec, 'Site', 'YearSR')
graphs <- lapply(nets, graph.incidence, weighted=TRUE)
save(graphs, nets, file="../data/nets.Rdata")


sp.lev <- calcSpec(nets, spec)
save(sp.lev, file='../data/sp.lev.Rdata')


