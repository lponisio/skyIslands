## setwd('~/Dropbox/skyislands')
rm(list=ls())
library(ggplot2, quietly = TRUE)
library(tidyverse, quietly = TRUE)
setwd("C:/")
source("lab_paths.R")
local.path

dir.bombus <- file.path(local.path, "skyIslands")
setwd(dir.bombus)
setwd('analysis/role')

source('src/initialize_beta.R')

## ************************************************************
## beta diversity as variation between sites,
## centroid for each year
## ************************************************************

dis <- mapply(function(a, b, c, d)
    calcBeta(comm= a, ## observed communities
             dis.method, ## dissimilarity metric
             nulls=b, ## null communities
             occ=binary, ## binary or abundance weighted?
             sub=type,
             zscore=FALSE), ## use Chase method not zscores
    a=comm$comm,
    b= nulls,
    SIMPLIFY=FALSE)

beta.dist <- makeBetaDataPretty()

save(beta.dist, file=sprintf("saved/results/pol_partnerVar_%s.Rdata",
                             net.type))

## *********************************************************
## Create summary figures
## *********************************************************

## Plants
## Beta calculated between years within a site (temporal)
load('saved/results/plant_partnerVar_Year.Rdata')

beta.dist %>% 
  ggplot(
  aes(y = reorder(GenusSpecies, dist), x = dist, color = Site)) +
  geom_point() +
  labs(x = "Partner variability", y = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(paste("Plant species", length(unique(beta.dist$GenusSpecies)))) #22 species
unique(beta.dist$GenusSpecies)

unique.species <- beta.dist %>%
  separate(GenusSpecies, into = c("Genus", "Species"), sep = " ") %>% 
  distinct(Genus, Species)

write.csv(unique.species, file = "saved/traits/plantSpecies_partner.csv", row.names = FALSE)

#Beta calculated between sites within a year (spatial)
load('saved/results/plant_partnerVar_Site.Rdata')

beta.dist %>% 
  ggplot(
    aes(y = reorder(GenusSpecies, dist), x = dist)) +
  geom_point() +
  labs(x = "Partner variability", y = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(paste("Plant species", length(unique(beta.dist$GenusSpecies)))) #48 species
unique(beta.dist$GenusSpecies)


## Pollinators
## Beta calculated between years within a site (temporal)
load('saved/results/pol_partnerVar_Year.Rdata')

beta.dist %>% 
  ggplot(
    aes(y = reorder(GenusSpecies, dist), x = dist)) +
  geom_point() +
  labs(x = "Partner variability", y = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(paste("Plant species", length(unique(beta.dist$GenusSpecies)))) #58 species
unique(beta.dist$GenusSpecies)

unique.species <- beta.dist %>%
  separate(GenusSpecies, into = c("Genus", "Species"), sep = " ") %>% 
  distinct(Genus, Species)

write.csv(unique.species, file = "saved/traits/polSpecies_partner.csv", row.names = FALSE)

#Beta calculated between sites within a year (spatial)
load('saved/results/pol_partnerVar_Site.Rdata')

beta.dist %>% 
  ggplot(
    aes(y = reorder(GenusSpecies, dist), x = dist)) +
  geom_point() +
  labs(x = "Partner variability", y = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(paste("Plant species", length(unique(beta.dist$GenusSpecies)))) #47 species
unique(beta.dist$GenusSpecies)

