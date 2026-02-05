##########################################
## Role/partner variability ~ plant pol traits
##########################################
rm(list=ls())
setwd("C:/")
source("lab_paths.R")
local.path
dir.bombus <- file.path(local.path, "skyIslands")
setwd(dir.bombus)
setwd('analysis/role')

library(tidyverse)
library(ggplot2)

## pollinators
load('saved/results/pol_partnerVar_Year.Rdata')

pol.traits <- read.csv('saved/traits/western_us_bee_traits.csv')
pols.SI <- read.csv('saved/traits/polSpecies_partner.csv')

colnames(pol.traits)
colnames(pols.SI)


pols.SI <- pols.SI %>% 
  mutate(GenusSpecies = paste(Genus, Species))


pol.traits <- pol.traits %>% 
  filter(GenusSpecies %in% pols.SI$GenusSpecies)

pol.beta.traits <- beta.dist %>% 
  left_join(pol.traits, by = "GenusSpecies")

pol.beta.traits[pol.beta.traits == ""] <- NA

write.csv(pol.beta.traits, file = "saved/traits/pol_beta_traits.csv")

## plants
load('saved/results/plant_partnerVar_Year.Rdata')

plant.traits <- read.csv('saved/traits/plantSpecies_partner.csv')

plant.traits <- plant.traits %>% 
  rename(GenusSpecies = genusSpecies)

plant.beta.traits <- beta.dist %>% 
  left_join(plant.traits, by = "GenusSpecies")

plant.beta.traits[plant.beta.traits == ""] <- NA

write.csv(plant.beta.traits, file = "saved/traits/plant_beta_traits.csv")

## plot
colnames(pol.beta.traits)
hist(pol.beta.traits$ForageDist_km)
hist(pol.beta.traits$MeanITD)

pol.beta.traits %>% 
  ggplot(aes(x = MeanITD, y = dist)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE)

pol.beta.traits %>% 
  ggplot(aes(x = ForageDist_km, y = dist)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE)

pol.beta.traits %>% 
  ggplot(aes(x = PollenCarry, y = dist)) +
  geom_boxplot()


colnames(plant.beta.traits)

plant.beta.traits %>% 
  ggplot(aes(x = symmetry, y = dist)) +
  geom_boxplot()

plant.beta.traits %>% 
  ggplot(aes(x = bloomPeriod, y = dist)) +
  geom_point() + 
  geom_smooth(method= "lm", se = TRUE)

