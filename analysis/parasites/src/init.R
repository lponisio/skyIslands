library(ggplot2)
library(brms)
library(dplyr)
library(tidyr)
library(tidybayes)
library(ggthemes)
library(car)
library(lme4)
library(bayestestR)


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

spec.net <- merge(spec.net, site.sum, all.x=TRUE)

traits <-
    read.csv("../../../skyIslands_saved/data/raw/bee_traits.csv")
traits$GenusSpecies <- fix.white.space(traits$GenusSpecies)
traits <- traits[, c("GenusSpecies", "Sociality", "Lecty", "MeanITD"),]

net.traits <- read.csv("../../data/networks_traits.csv")
net.traits <- net.traits[, c("GenusSpecies", "r.degree"),]

traits <- merge(traits, net.traits, by="GenusSpecies", all.x=TRUE)

spec.net <- merge(spec.net, traits, all.x=TRUE, by="GenusSpecies")

dir.create(path="saved", showWarnings = FALSE)
dir.create(path="saved/tables", showWarnings = FALSE)

spec.net <- spec.net[order(spec.net$Site),]

## create a dumby varaible "Weight" to deal with the data sets being at
## different levels to get around the issue of having to pass in one
## data set into brms
spec.net$YearSR <- paste(spec.net$Year, spec.net$SampleRound, sep=";")

## will need to modify when we have multiple years
spec.net <- makeDataMultiLevel(spec.net, "Site", "YearSR")

spec.net[, variables.to.log] <- log(spec.net[,variables.to.log] + 1)

##  center all of the x variables, need to use unique values to avoid
##  repetition by the number of specimens
unique.site.vals <- spec.net[spec.net$Weights == 1, ]
unique.site.vals[, vars] <- apply(unique.site.vals[, vars], 2, standardize)
spec.net[, vars] <- NULL

unique.site.vals <- unique.site.vals[, c("Site", "YearSR", vars)]

dim(spec.net)
spec.net <- merge(spec.net, unique.site.vals)
dim(spec.net)

## create a dumby varaible "WeightPar" for the parasite data. The
## original intention was to keep stan from dropping data for
## site-level models, but weight is 0 for parasite models.

spec.net$WeightsPar <- 1
spec.net$WeightsPar[is.na(spec.net$ParasitePresence) | spec.net$Apidae != 1] <- 0

## stan drops all NA data, so can set ParasitePresence to 0 with WeightsPar
## to keep it in the models
spec.net$ParasitePresence[is.na(spec.net$ParasitePresence | spec.net$Apidae != 1)] <- 0

spec.all <- spec.net
spec.all[, vars] <- apply(spec.all[, vars], 2, standardize)

## bombus only data
spec.bombus <- spec.all
spec.bombus$WeightsPar[spec.bombus$Genus != "Bombus"] <- 0

## apis only data
spec.apis <- spec.all
spec.apis$WeightsPar[spec.apis$Genus != "Apis"] <- 0

## melissodes only data
spec.melissodes <- spec.all
spec.melissodes$WeightsPar[spec.melissodes$Genus != "Melissodes"] <- 0

spec.apidae <- spec.all
spec.apidae$WeightsPar[spec.melissodes$Family != "Apidae"] <- 0

## not enough replication of species in Halictidae
spec.all$lWeightsPar[spec.all$Family =="Halictidae"] <- 0
## not enough replication of species in Colletidae
spec.all$lWeightsPar[spec.all$Family =="Colletidae"] <- 0
