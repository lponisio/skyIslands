rm(list=ls())
library(igraph)
library(vegan)
library(fields)
library(bipartite)
setwd("~/Dropbox/skyIslands/dataPrep")
source("src/misc.R")
source("src/prepNets.R")
source("src/specialization.R")

src.dir <- '../../skyIslands_saved/data/relational/relational/traditional/'
spec <-
    read.csv(file.path(src.dir, "specimens-complete.csv"))

## did not complete full sampling rounds in any of these sites. Was
## just scouting.
site.2.drop <- c("JM", "SS", "UK", "CC")
spec <- spec[!spec$Site %in% site.2.drop,]
spec <- droplevels(spec)

## get specimen data ready
spec$GenusSpecies <- fix.white.space(paste(spec$Genus,
                          spec$Species,
                          spec$SubSpecies))

spec$PlantGenusSpecies <-  fix.white.space(paste(spec$PlantGenus,
                                          spec$PlantSpecies,
                                          spec$PlantVar, spec$PlantSubSpecies))

spec$Int <-  fix.white.space(paste(spec$GenusSpecies, spec$PlantGenusSpecies))
spec$IntGen <-  fix.white.space(paste(spec$Genus, spec$PlantGenus))

spec <- spec[spec$PlantGenusSpecies != "",]
spec <- spec[spec$GenusSpecies != "",]

## drop non-bees
spec <- spec[spec$Family %in% c("Andrenidae", "Apidae",
                                "Colletidae", "Halictidae", "Megachilidae"),]

spec$Date <- as.Date(spec$Date, format='%m/%d/%y')
spec$Doy <- as.numeric(strftime(spec$Date, format='%j'))
spec$Year <- as.numeric(format(spec$Date,'%Y'))


## drop the the 2017 sample of PL because it was on fire for other
## sampling rounds
spec <- spec[!(spec$Site == "PL" & spec$Year == "2017"),]


## calculate network metric of interact
## make a plant by pol matrix
prep.comm <- aggregate(spec$GenusSpecies,
                       list(site=spec$PlantGenusSpecies,
                            sp=spec$GenusSpecies), length)
comm <- samp2site.spp(prep.comm$site, prep.comm$sp, prep.comm$x,
                      FUN=mean)

## calculate specialization
d <- specieslevel(comm, index="d")

spec$PolSpec <- d$'higher level'$'d'[match(spec$GenusSpecies,
                                             rownames(d$'higher level'))]

spec$PlantSpec <- d$'lower level'$'d'[match(spec$PlantGenusSpecies,
                                             rownames(d$'lower level'))]

save(spec, file="../data/spec.Rdata")


### networks
nets <- breakNet(spec, 'Site', 'Year')
graphs <- lapply(nets, graph.incidence, weighted=TRUE)
save(graphs, nets, file="../data/nets.Rdata")

sp.lev <- calcSpec(nets, spec)
save(sp.lev, file='../data/sp.lev.Rdata')


