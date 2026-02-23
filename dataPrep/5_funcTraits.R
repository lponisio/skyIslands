rm(list=ls())
library(FD)
library(vegan)
library(readxl)

setwd("~/")
setwd('C:/')
source("lab_paths.R")
local.path

library(tidyverse)
dir.bombus <- file.path(local.path, "skyIslands/dataPrep")
setwd(dir.bombus)

## calculate trait uniqueness and originality based on Coux et
## al. 2016 for each site-year-round

source("src/calcFuncUniqOrig.R")
source("src/misc.R")
load('../data/spec_net.Rdata')
## only working with bees
bee.families <- c("Andrenidae", "Apidae", "Colletidae", "Halictidae",
                  "Megachilidae")

spec.net <- spec.net[spec.net$Family %in% bee.families,]

traits <-
  read_excel("../../PonisioLab/bee_traits/western_us_bee_traits.xlsx",
             sheet=1)
traits <- traits[traits$GenusSpecies %in% spec.net$GenusSpecies,]
traits <- as.data.frame(traits)

missing <- spec.net[!spec.net$GenusSpecies %in%
                                 traits$GenusSpecies,
                    c("GenusSpecies", "Year")]
missing

itd <- traits %>%
  group_by(Genus) %>%
  summarize(MeanITD = mean(MeanITD, na.rm=TRUE))

traits$MeanITD[is.na(traits$MeanITD)] <- itd$MeanITD[
  match(traits$Genus[is.na(traits$MeanITD)], itd$Genus)]

bee.traits <- c("NestLocation", "NestPartitions","NestConstruction",
                "NestLinning", "PrimaryNestMaterial", "MeanITD",
                "PollenCarry", "ReproStrat", "Sociality"
                )

traits <- traits[, c("GenusSpecies", bee.traits)]

traits[traits == ""] <- NA


rownames(traits) <- traits$GenusSpecies
traits$GenusSpecies <- NULL

## We grouped the traits based on what they discuss. so for exaple, we
## have 3 nesting traits, so weighted each one by 1/3.

bee.weights <- c(rep(1/5, 5), rep(1, 2), rep(1/2, 2))
length(bee.weights) == length(bee.traits)


# make the three IDs once
spec.net <- spec.net %>%
  tidyr::unite("SiteYearSR", Site, Year, SampleRound, sep = ";", remove = FALSE) %>%
  tidyr::unite("SiteYear",   Site, Year,           sep = ";", remove = FALSE)

# run all three with one pattern
spec.net <- add_func_uniq_orig(spec.net, traits,
                               id_col = "Site",
                               traits.2.keep = bee.traits, weights = bee.weights,
                               suffix = "Site",
                               add_fd = TRUE)

spec.net <- add_func_uniq_orig(spec.net, traits,
                               id_col = "SiteYearSR",
                               traits.2.keep = bee.traits, weights = bee.weights,
                               suffix = "SiteYearSR",
                               add_fd = TRUE)

spec.net <- add_func_uniq_orig(spec.net, traits,
                               id_col = "SiteYear",
                               traits.2.keep = bee.traits, weights = bee.weights,
                               suffix = "SiteYear",
                               add_fd = TRUE)

save(spec.net, file="../data/spec_traits.Rdata")
write.csv(spec.net, file="../data/spec_traits.csv", row.names=FALSE)

# Merges network traits by site and round to spec.net
load('../data/sp_network_mets_sr_pretty.RData')

sp.network.metrics <- sp.network.metrics %>%
  select(GenusSpecies, Site, Year, SampleRound, zdegree,
         zweighted.betweenness, zweighted.closeness, zd,
         normalised.degree)

traits$GenusSpecies <- rownames(traits)
rownames(traits) <- NULL

dim(spec.net)
spec.net <- merge(spec.net, sp.network.metrics, all.x=TRUE)
dim(spec.net)
spec.net <- merge(spec.net, traits, all.x=TRUE)
dim(spec.net)

save(spec.net, file="../data/spec_traits.Rdata")
write.csv(spec.net, file="../data/spec_traits.csv", row.names=FALSE)

# ## Merges network traits by year to spec.net
# load('../data/sp_network_mets_Year_pretty.RData')
#  
# sp.network.metrics <- sp.network.metrics %>%
#    select(GenusSpecies, Site, Year, SampleRound, zdegree,
#           zweighted.betweenness, zweighted.closeness, zd,
#           normalised.degree)
#  
# traits$GenusSpecies <- rownames(traits)
# rownames(traits) <- NULL
#  
# dim(spec.net)
# spec.net <- merge(spec.net, sp.network.metrics, all.x=TRUE)
# dim(spec.net)
# spec.net <- merge(spec.net, traits, all.x=TRUE)
# dim(spec.net)
#  
# save(spec.net, file="../data/spec_traits_year.Rdata")
# write.csv(spec.net, file="../data/spec_traits_year.csv", row.names=FALSE)

spec.bee.microbes <- spec.net %>%
  filter(Family != "Syrphidae") %>%
  select(UniqueID, TempID, GenusSpecies, Site, Year, SampleRound,
         PlantGenusSpecies, NestLocation, NestPartitions,
         NestConstruction, PrimaryNestMaterial, MeanITD, PollenCarry,
         ReproStrat, Sociality,
         originalitySiteYearSR, uniqSiteYearSR,
         originalitySiteYear, uniqSiteYear,
         originalitySite, uniqSite,
         zdegree,
         zweighted.betweenness, zweighted.closeness, zd,
         normalised.degree)


write.csv(spec.bee.microbes,
          file="../data/si_bee_microbes_traits.csv", row.names=FALSE)
