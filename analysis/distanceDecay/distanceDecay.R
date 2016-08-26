rm(list=ls())
library(vegan)
library(fields)
setwd("~/Dropbox/SkyIslands/analysis/distanceDecay")
source('src/misc.R')
source('src/distDecay.R')

geo <-
  read.csv("../../data/relational/data/relational/tables/geography.csv")

spec <-
  read.csv("../data/spec.csv")

##distance dissimilarity
dist.site <- rdist.earth(cbind(geo$Long, geo$Lat),
                         cbind(geo$Long, geo$Lat))

c.dist <- dist.site[lower.tri(dist.site)]

## dissimilarity of plants, pol, int

sp <- ppint.dis(spec, c.dist=c.dist)

sp.raup <- ppint.dis(spec, c.dist=c.dist,
                     method.dist="raup")

gen <- ppint.dis(spec, types=
                 c("Genus", "PlantGenus", "IntGen"),
                 c.dist=c.dist,
                 sub="all")

## just bees
spec.bee <- spec[spec$Order == "Hymenoptera",]
sp.bee <- ppint.dis(spec.bee, types= c("GenSp", "PlantGenSp", "Int"),
                    c.dist=c.dist, sub="bee")
gen.bee <- ppint.dis(spec.bee, types= c("Genus", "PlantGenus", "IntGen"),
                     c.dist=c.dist, sub="bee")

## just leps
spec.lep <- spec[spec$Order == "Lepidoptera",]
sp.lep <- ppint.dis(spec.lep, types= c("GenSp", "PlantGenSp", "Int"),
                    c.dist=c.dist, sub="lep")
gen.lep <- ppint.dis(spec.lep, types= c("Genus", "PlantGenus", "IntGen"),
                     c.dist=c.dist, sub="lep")

## just syrphids
spec.syr <- spec[spec$Family == "Syrphidae",]
sp.syr <- ppint.dis(spec.syr, types= c("GenSp", "PlantGenSp", "Int"),
                    c.dist=c.dist, sub="syrphid")
gen.syr <- ppint.dis(spec.syr, types= c("Genus", "PlantGenus", "IntGen"),
                     c.dist=c.dist, sub="syrphid")
