## setwd('~/Dropbox/skyislands')
rm(list=ls())
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

save(beta.dist, file=sprintf("saved/results/partnerVar_%s.Rdata",
                             net.type))

## *********************************************************
## Create summary figures
## *********************************************************

library(ggplot2)
library(tidyverse)

#Beta calculated between years within a site (temporal)
#source('saved/results/partnerVar_Year.Rdata')

unique(beta.dist$Site)

beta.dist %>% 
  ggplot(aes(x = dist)) +
  geom_histogram() +
  facet_wrap(Site~Year)

beta.dist %>% 
  filter(Site == 'CH') %>% 
  ggplot(
    aes(x = dist)) +
  geom_histogram() +
  facet_wrap(~ Year)


beta.dist %>% 
  ggplot(
  aes(y = fct_reorder(GenusSpecies, dist), x = dist)) +
  geom_point() +
  labs(x = "Partner variability", y = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(paste("Plant species", length(unique(beta.dist$GenusSpecies)))) #22 species
unique(beta.dist$GenusSpecies)

#Beta calculated between sites within a year (spatial)
#source('saved/results/partnerVar_Site.Rdata')

unique(beta.dist$Year)

beta.dist %>% 
  ggplot(
    aes(x = dist)) +
  geom_histogram() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(Year ~ Site)

beta.dist %>% 
  filter(Year == 2018) %>% 
  ggplot(
    aes(x = dist)) +
  geom_histogram() +
  facet_wrap( ~ Site)



print(paste("Plant species", length(unique(beta.dist$GenusSpecies)))) #48 species
unique(beta.dist$GenusSpecies)

