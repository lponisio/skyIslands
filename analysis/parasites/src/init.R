library(ggplot2)
library(brms)
library(dplyr)
library(tidyr)
library(tidybayes)
library(ggthemes)
library(bayestestR)

dir.create(path="saved", showWarnings = FALSE)
dir.create(path="saved/tables", showWarnings = FALSE)

load('../../data/spec_net.Rdata')
site.sum <- read.csv("../../data/sitestats.csv")

spec.net <- spec.net[!is.na(spec.net$GenusSpecies),]

spec.net <- spec.net[spec.net$Family != "Syrphidae",]

parasites <- c(#"AspergillusSpp", ## problematic parasite!
  "AscosphaeraSpp",
  "ApicystisSpp",
  "CrithidiaExpoeki",
  "CrithidiaBombi",
  "CrithidiaSpp",
  "NosemaBombi",
  "NosemaCeranae")

## Merging specimen data with site characteristic data.
print("Before merge with site characteristics")
print(dim(spec.net))
spec.net <- merge(spec.net, site.sum, all.x=TRUE)
print("After merge with site characteristics")
print(dim(spec.net))

## Merging specimen data with trait data and network trait data.
traits <-
    read.csv("../../../skyIslands_saved/data/raw/bee_traits.csv")
traits$GenusSpecies <- fix.white.space(traits$GenusSpecies)
traits <- traits[, c("GenusSpecies", "Sociality", "Lecty", "MeanITD"),]

net.traits <- read.csv("../../data/networks_traits.csv")
net.traits <- net.traits[, c("GenusSpecies", "r.degree"),]
print("Before merge with network traits to species traits")
print(dim(traits))
traits <- merge(traits, net.traits, by="GenusSpecies", all.x=TRUE)
print("After merge with network traits to species traits")
print(dim(traits))
## Merge all the trait data to the specimen data.
spec.net <- merge(spec.net, traits, all.x=TRUE, by="GenusSpecies")
print("After merge with specimen traits")
print(dim(spec.net))
spec.orig <- spec.net

## Make SEM weights and standardize data.
spec.net <- prepDataSEM(spec.net, variables.to.log, variables.to.log.1, 
                        vars_yearsr, vars_sp)

## bombus only data
spec.bombus <- spec.net
spec.bombus$WeightsPar[spec.bombus$Genus != "Bombus"] <- 0

## apis only data
spec.apis <- spec.net
spec.apis$WeightsPar[spec.apis$Genus != "Apis"] <- 0

## melissodes only data
spec.melissodes <- spec.net
spec.melissodes$WeightsPar[spec.melissodes$Genus != "Melissodes"] <- 0

spec.apidae <- spec.net
spec.apidae$WeightsPar[spec.apidae$Family != "Apidae"] <- 0


