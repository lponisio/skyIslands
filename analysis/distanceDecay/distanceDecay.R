## setwd("~/Dropbox/skyIslands/")
rm(list=ls())
setwd("analysis/distanceDecay")
source("src/initialize.R")

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
