rm(list=ls())
library(igraph)
library(vegan)
library(fields)
library(bipartite)
setwd("~/Dropbox/skyIslands/analysis/data")
source("src/misc.R")
source("src/prepNets.R")
source("src/specialization.R")
src.dir <- '../../data/relational/data/relational/traditional/'
spec <-
  read.csv(file.path(src.dir, "specimens-complete.csv"))

## get specimen data ready
spec <-  spec[spec$Species != "",]
spec$GenusSpecies <- fix.white.space(paste(spec$Genus, 
                          spec$Species,
                          spec$SubSpecies))

spec$PlantGenusSpecies <-  fix.white.space(paste(spec$PlantGenus,
                                          spec$PlantSpecies,
                                          spec$PlantVar, spec$PlantSubSpecies))
spec$Int <-  fix.white.space(paste(spec$GenusSpecies, spec$PlantGenusSpecies))
spec$IntGen <-  fix.white.space(paste(spec$Genus, spec$PlantGenus))

prep.comm <- aggregate(spec$GenusSpecies,
                       list(site=spec$PlantGenusSpecies,
                            sp=spec$GenusSpecies), length)
comm <- samp2site.spp(prep.comm$site, prep.comm$sp, prep.comm$x)

d <- specieslevel(comm, index="d")

spec$PolSpec <- d$'higher level'$'d'[match(spec$GenusSpecies,
                                             rownames(d$'higher level'))]

spec$PlantSpec <- d$'lower level'$'d'[match(spec$PlantGenusSpecies,
                                             rownames(d$'lower level'))]

spec$Date <- as.Date(spec$Date, format='%m/%d/%y')
spec$Doy <- as.numeric(strftime(spec$Date, format='%j'))
spec$Year <- as.numeric(format(spec$Date,'%Y'))
                                         
write.csv(spec, "spec/spec.csv", row.names=FALSE)
save(spec, file="spec/spec.Rdata")


### networks
nets <- breakNet(spec, 'Site', 'Year')
graphs <- lapply(nets, graph.incidence, weighted=TRUE, add.names=NA)
save(graphs, nets, file="networks/nets.Rdata")

sp.lev <- calcSpec(nets, spec)
save(sp.lev, file='spec/sp.lev.Rdata')


