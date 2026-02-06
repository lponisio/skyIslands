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
source('src/initialize_nulls.R')

library(tidyverse, quietly = TRUE)
library(ggplot2, quietly = TRUE)

## Calculate bloom and flight period
phen <- spec %>%
  filter(Sex == 'f') %>% 
  mutate(PlantGenusSpecies = paste(PlantGenus, PlantSpecies)) %>% 
  group_by(PlantGenusSpecies, Year) %>%
  mutate(
    bloom_start = min(Doy),
    bloom_end = max(Doy),
    bloom_period = bloom_end - bloom_start)

phen <- phen %>% 
  group_by(GenusSpecies, Year) %>% 
  mutate(
    beeEmergence_start = min(Doy),
    beeEmergence_end = max(Doy),
    flight_period = beeEmergence_end - beeEmergence_start)

#average across all columns
phen.av <- phen %>% 
  group_by(GenusSpecies, Site, Year) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop") %>% 
  mutate(GenusSpeciesSiteYear = paste0(GenusSpecies, Site, Year))


#######################################
## pollinators
#######################################
load('saved/results/pol_partnerVar_Year.Rdata')

pol.traits <- read.csv('saved/traits/western_us_bee_traits.csv')
pols.SI <- read.csv('saved/traits/polSpecies_partner.csv')

pols.SI <- pols.SI %>% 
  mutate(GenusSpecies = paste(Genus, Species))


pol.traits <- pol.traits %>% 
  filter(GenusSpecies %in% pols.SI$GenusSpecies)

pol.beta.traits <- beta.dist %>% 
  left_join(pol.traits, by = "GenusSpecies") %>% 
  mutate(GenusSpeciesSiteYear = paste0(GenusSpecies, Site, Year))

pol.beta.traits[pol.beta.traits == ""] <- NA

#right join (keeps all rows from pol.beta.traits)
pol.beta.traits.phen <- phen.av %>% 
  right_join(pol.beta.traits, by = "GenusSpeciesSiteYear")

dim(pol.beta.traits)
dim(phen)
dim(pol.beta.traits.phen)

write.csv(pol.beta.traits.phen, file = "saved/traits/pol_beta_traits.csv")

colnames(pol.beta.traits.phen)


pol.beta.traits.phen %>% 
  ggplot(aes(x = nestedrank, y = dist)) +
  geom_point() +
  geom_smooth(method = "lm")








######################################
## plants
######################################
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

