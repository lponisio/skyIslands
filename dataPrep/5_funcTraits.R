rm(list=ls())
library(FD)
library(vegan)

setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyIslands')
setwd('dataPrep')

#calculate trait uniqueness and originality based on Coux et al. 2016
source("src/calcFuncUniqOrig.R")
source("src/misc.R")

load('../data/spec.Rdata')

traits <- read.csv("../../skyIslands_saved/data/raw/bee_traits.csv")
traits$GenusSpecies <- fix.white.space(traits$GenusSpecies)

traits <- traits[traits$GenusSpecies %in% spec$GenusSpecies,]

bee.traits <- c("NestLocation", "NestMaterial","NestConstruction","Lecty",
                 "MeanITD")

rownames(traits) <- traits$GenusSpecies
traits$GenusSpecies <- NULL

## We grouped the traits based on what they are related to. For
## example, we have 3 nesting-related traits, so weighted each one by
## 1/3. This will avoid overweighting any trait type.

bee.weights <- c(rep(1/2, 3), rep(1, 2))

## *****************************************************************
## species-level function diversity metrics
## *****************************************************************
## calculates data for two new columns we are adding to database,
## "originality" and "uniq"
bee.func <- calcFuncUniqOrig(traits,
                             traits.2.keep=bee.traits,
                             weights=bee.weights,
                             type="all spec")

bee.func$GenusSpecies <- rownames(bee.func)
traits$GenusSpecies <- rownames(traits)
rownames(traits) <- NULL
rownames(bee.func) <- NULL

traits <- merge(traits, bee.func, by="GenusSpecies")

all.traits <- colnames(traits)[!colnames(traits) %in%
                               c("GenusSpecies")]
## save prepped data
spec <- cbind(spec, traits[, all.traits][match(spec$GenusSpecies,
                           traits$GenusSpecies),])

## spec$r.degree <- as.numeric(spec$r.degree)

print(paste("missing trait data",
            sort(unique(spec$GenusSpecies[is.na(spec$Lecty)]))))

write.csv(traits, file="../data/traits.csv",
          row.names=FALSE)

save(spec, file="../data/spec_traits.Rdata")
write.csv(spec, file="../data/spec_traits.csv", row.names=FALSE)

## *****************************************************************
## site-level function diversity metrics
## *****************************************************************
