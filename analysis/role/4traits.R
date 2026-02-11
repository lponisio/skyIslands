#####################################################
## -- This script merges trait data with        -- ##
## -- role/partner variability and climate data -- ##
#####################################################
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


pol.traits <- read.csv(file = 'saved/traits/pol_beta_traits.csv')
plant.traits <- read.csv(file = 'saved/traits/plant_beta_traits.csv')
load('saved/results/pol_partnerVar_Site.Rdata')

dim(spec.net)
colnames(pol.traits)


## --- calculate emergence and flight period -- ##
## -- average across columns species within a site and year --- ##
spec <- spec.net %>% 
  group_by(GenusSpecies, Year) %>% 
  mutate(beeEmergence_start = min(Doy),
  beeEmergence_end = max(Doy),
  flight_period = beeEmergence_end - beeEmergence_start,
    .groups = "drop") %>% 
  group_by(GenusSpecies, Site, Year) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), 
            .groups = "drop") %>% 
  mutate(GenusSpeciesSiteYear = paste0(GenusSpecies, Site, Year))

pol.beta <- beta.dist %>%
  mutate(GenusSpeciesSiteYear = paste0(GenusSpecies, Site, Year))

## --- merge species level network metrics, role variability, trait, and climate data --- ##
pol.beta.traits <- spec %>%
  right_join(pol.beta, by = "GenusSpeciesSiteYear") # keeps all rows from pol.beta

write.csv(pol.beta.traits, file = 'saved/traits/pol_beta_traits.csv')

## --- merge climate data --- ##
load('saved/traits/climateVariability.csv')

colnames(pol.traits)

pol.trait.climate <- pol.traits %>%
  left_join(climateVar, by = c("Site", "Year"))

write.csv(pol.trait.climate, file = 'saved/traits/pol_trait_climate.csv')


## -- Assess sampling start date for each site by year -- ##
spec.net %>% 
  ggplot(aes(x = Doy, y = SampleRound)) +
  geom_point() +
  facet_wrap(Site~Year)


########################################################################
## -- Merge plant traits, network metric, and role variability data-- ##
########################################################################
load('saved/results/plant_partnerVar_Year.Rdata')

## Calculate bloom and flight period
phen <- spec %>%
  filter(Sex == 'f') %>%
   mutate(PlantGenusSpecies = paste(PlantGenus, PlantSpecies)) %>%
   group_by(PlantGenusSpecies, Year) %>%
   mutate(
     bloom_start = min(Doy),
     bloom_end = max(Doy),
     bloom_period = bloom_end - bloom_start)

plant.traits <- read.csv('saved/traits/plantSpecies_partner.csv')

plant.traits <- plant.traits %>%
   rename(GenusSpecies = genusSpecies)

plant.beta.traits <- beta.dist %>%
   left_join(plant.traits, by = "GenusSpecies")

plant.beta.traits[plant.beta.traits == ""] <- NA

write.csv(plant.beta.traits, file = "saved/traits/plant_beta_traits.csv")

## -- Merge climate variability data -- ##





# #######################################
# ## -- pollinators trait data -- ##
# #######################################
# load('saved/results/pol_partnerVar_Year.Rdata')
# 
# pol.traits <- read.csv('saved/traits/western_us_bee_traits.csv')
# pols.SI <- read.csv('saved/traits/polSpecies_partner.csv')
# 
# pols.SI <- pols.SI %>% 
#   mutate(GenusSpecies = paste(Genus, Species))
# 
# 
# pol.traits <- pol.traits %>% 
#   filter(GenusSpecies %in% pols.SI$GenusSpecies)
# 
# pol.beta.traits <- beta.dist %>% 
#   left_join(pol.traits, by = "GenusSpecies") %>% 
#   mutate(GenusSpeciesSiteYear = paste0(GenusSpecies, Site, Year))
# 
# pol.beta.traits[pol.beta.traits == ""] <- NA
# 
# #right join (keeps all rows from pol.beta.traits)
# pol.beta.traits.phen <- phen.av %>% 
#   right_join(pol.beta.traits, by = "GenusSpeciesSiteYear")
# 
# dim(pol.beta.traits)
# dim(phen)
# dim(pol.beta.traits.phen)
# 
# write.csv(pol.beta.traits.phen, file = "saved/traits/pol_beta_traits.csv")
# 



