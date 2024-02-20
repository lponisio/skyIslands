m(list=ls())
library(FD)
library(vegan)

source("lab_paths.R")
local.path

dir.bombus <- file.path(local.path, "skyIslands/dataPrep")
setwd(dir.bombus)

## calculate trait uniqueness and originality based on Coux et al. 2016
## currently traits are across all SI, need to update to be meadow
## specific

## also many bees are missing traits, need to fill in

source("src/calcFuncUniqOrig.R")
source("src/misc.R")
load('../data/spec_net.Rdata')

traits <- read.csv("../../skyIslands_saved/data/raw/bee_traits.csv")
traits$GenusSpecies <- fix.white.space(traits$GenusSpecies)

traits <- traits[traits$GenusSpecies %in% spec.net$GenusSpecies,]

bee.traits <- c("NestLocation",
                "PrimaryNestMaterial","NestConstruction",
                "Lecty",
                 "MeanITD", "Sociality")

m.traits<- traits$GenusSpecies[apply(traits[, bee.traits], 1, 
                           function(x) any(is.na(x))
                           )]

write.csv(m.traits,
          file="../../skyIslands_saved/data/checks/missing_traits.csv")

rownames(traits) <- traits$GenusSpecies
traits$GenusSpecies <- NULL


## We grouped the traits based on what they are related to. For
## example, we have 3 nesting-related traits, so weighted each one by
## 1/3. This will avoid overweighting any trait type.

bee.weights <- c(rep(1/3, 3), rep(1, 3))

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
spec.net <- cbind(spec.net, traits[, all.traits][match(spec.net$GenusSpecies,
                           traits$GenusSpecies),])

## spec.net$r.degree <- as.numeric(spec.net$r.degree)

print(paste("missing trait data",
            sort(unique(spec.net$GenusSpecies[is.na(spec.net$Lecty)]))))

write.csv(traits, file="../data/traits.csv",
          row.names=FALSE)

save(spec.net, file="../data/spec_traits.Rdata")
write.csv(spec.net, file="../data/spec_traits.csv", row.names=FALSE)

## *****************************************************************
## site-level function diversity metrics
## *****************************************************************
