library(ggplot2)
library(brms)
library(dplyr)
library(tidyr)
library(tidybayes)
library(ggthemes)
library(bayestestR)
library(car)

dir.create(path="saved", showWarnings = FALSE)
dir.create(path="saved/tables", showWarnings = FALSE)

load('../../data/spec_net.Rdata')

parasites <- c(
    "AscosphaeraSpp",
    "ApicystisSpp",
    "CrithidiaExpoeki",
    "CrithidiaBombi",
    "CrithidiaSpp",
    "CrithidiaMellificae",
    "NosemaBombi",
    "NosemaCeranae")

## ## Merging specimen data with site characteristic data.
## print("Before merge with site characteristics")
## print(dim(spec.net))
## spec.net <- merge(spec.net, site.sum, all.xy=TRUE)
## print("After merge with site characteristics")
## print(dim(spec.net))

## ## Merging species data with individual level data. 
## print("Before merge with ind level data")
## print(dim(spec.net))
## spec.net <- merge(spec.net, sp.sum, all.x= TRUE)
## print("After merge with ind level data")
## print(dim(spec.net))


## ## Merging specimen data with trait data and network trait data.
## traits <-
##     read.csv("../../../skyIslands_saved/data/raw/bee_traits.csv")
## traits$GenusSpecies <- fix.white.space(traits$GenusSpecies)
## traits <- traits[, c("GenusSpecies", "Sociality", "Lecty", "MeanITD"),]

## net.traits <- read.csv("../../data/networks_traits.csv")
## net.traits <- net.traits[, c("GenusSpecies", "r.degree"),]
## print("Before merge with network traits to species traits")
## print(dim(traits))
## traits <- merge(traits, net.traits, by="GenusSpecies", all.x=TRUE)
## print("After merge with network traits to species traits")
## print(dim(traits))
## ## Merge all the trait data to the specimen data.
## spec.net <- merge(spec.net, traits, all.x=TRUE, by="GenusSpecies")
## print("After merge with specimen traits")
## print(dim(spec.net))

to.drop <- !grepl("Pan_", colnames(spec.net))

spec.net <- spec.net[, to.drop]




