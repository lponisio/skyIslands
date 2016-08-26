rm(list=ls())
library(vegan)
library(fields)
library(bipartite)
setwd("~/Dropbox/SkyIslands/analysis/data")
source('src/samp2site_spp.R')
source("src/misc.R")
spec <-
  read.csv("../../data/relational/data/relational/traditional/specimens-complete.csv")

##get specimen data ready
spec <-  spec[spec$Species != "",]
spec$GenSp <- fix.white.space(paste(spec$Genus, 
                          spec$Species,
                          spec$SubSpecies))

spec$PlantGenSp <-  fix.white.space(paste(spec$PlantGenus,
                                          spec$PlantSpecies,
                                          spec$PlantVar, spec$PlantSubSpecies))
spec$Int <-  fix.white.space(paste(spec$GenSp, spec$PlantGenSp))
spec$IntGen <-  fix.white.space(paste(spec$Genus, spec$PlantGenus))

prep.comm <- aggregate(spec$GenSp,
                       list(site=spec$PlantGenSp,
                            sp=spec$GenSp), length)
comm <- samp2site.spp(prep.comm$site, prep.comm$sp, prep.comm$x)

d <- specieslevel(comm, index="d")

spec$PolSpec <- d$'higher level'$'d'[match(spec$GenSp,
                                             rownames(d$'higher level'))]

spec$PlantSpec <- d$'lower level'$'d'[match(spec$PlantGenSp,
                                             rownames(d$'lower level'))]
                                         
write.csv(spec, "spec.csv", row.names=FALSE)
save(spec, file="spec.Rdata")
