rm(list=ls())
library(FD)
library(vegan)
setwd("~/")
source("lab_paths.R")
local.path

library(tidyverse)
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
traits$Genus <- sapply(strsplit(traits$GenusSpecies, "\\ "),
                       function(x) x[1])

traits <- traits[traits$GenusSpecies %in% spec.net$GenusSpecies,]

itd <- traits %>%
  group_by(Genus) %>%
  summarize(MeanITD = mean(MeanITD, na.rm=TRUE))

traits$MeanITD[is.na(traits$MeanITD)] <- itd$MeanITD[
  match(traits$Genus[is.na(traits$MeanITD)], itd$Genus)]

bee.traits <- c("NestLocation", "PrimaryNestMaterial",
                "NestConstruction","NestPartitions","Sociality",
                "MeanITD")

m.traits <- traits$GenusSpecies[apply(traits[, bee.traits], 1, 
                           function(x) any(is.na(x))
                           )]

write.csv(m.traits,
          file="../../skyIslands_saved/data/checks/missing_traits.csv")


## an assumption to try to get things to run
traits$Sociality[traits$Genus == "Lasioglossum"] <- "solitary"

rownames(traits) <- traits$GenusSpecies
traits$GenusSpecies <- NULL

## check on unique values 
apply(traits[, bee.traits], 2, unique)

## We grouped the traits based on what they are related to. For
## example, we have 3 nesting-related traits, so weighted each one by
## 1/3. This will avoid overweighting any trait type.

bee.weights <- c(rep(1/4, 4), rep(1, 2))

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

spec.net$SiteYearSr <- paste(spec.net$Site, spec.net$Year,
                             spec.net$SampleRound,
                             sep=";")
## only working with bees
bee.families <- c("Andrenidae", "Apidae", "Colletidae", "Halictidae",
                  "Megachilidae")
spec.net <- spec.net[spec.net$Family %in% bee.families,]

spec.comm <- spec.net[, c("GenusSpecies", "SiteYearSr")]
spec.comm <- spec.comm[spec.comm$GenusSpecies != "",]

comms <- spec.comm %>% pivot_longer(cols = -SiteYearSr) %>%
  group_by(SiteYearSr,value) %>% summarise(N=n()) %>%
  pivot_wider(names_from = value, values_from=N) %>%
  replace(is.na(.),0)

comms <- as.data.frame(comms)
rownames(comms) <- comms$SiteYearSr
comms$SiteYearSr <- NULL

colnames(comms)[!colnames(comms) %in% rownames(traits)]

bee.func <- calcFuncUniqOrig(traits,
                             traits.2.keep=bee.traits,
                             weights=bee.weights,
                             type="networks",
                             a = comms,
                             w.abun=TRUE
)
