####################################################################
## Calculates phenology for each plant species within a single year
####################################################################
rm(list=ls())
setwd("C:/")
source("lab_paths.R")
local.path

dir.bombus <- file.path(local.path, "skyIslands")
setwd(dir.bombus)
setwd('analysis/role')
source('src/initialize_nulls.R')
library(tidyverse, quitely = TRUE)

phen <- spec %>%
  mutate(PlantGenusSpecies = paste(PlantGenus, PlantSpecies)) %>% 
  group_by(PlantGenusSpecies, Year) %>%
  mutate(
    bloom_start = min(Doy),
    bloom_end = max(Doy))

avg.phen <- phen %>% 
  group_by(PlantGenusSpecies, Year) %>% 
  summarize(
    av.bloom.start = mean(bloom_start),
    av.bloom.end = mean(bloom_end),
    n = n(),
    .groups = "drop")
