## setwd('~/Dropbox/skyislands')
rm(list=ls())
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

#Beta dis calculated between years within a site
#source('saved/results/partnerVar_Year.Rdata')

beta.dist %>% 
  filter(Site == '') +
  ggplot(
    aes(x = dist)) +
      geom_histogram() +
  facet_wrap(~ Year)


beta.dist %>% 
  ggplot(
  aes(x = GenusSpecies, y = dist)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Year)


print(paste("Plant species", length(unique(beta.dist$GenusSpecies)))) #22 species
unique(beta.dist$GenusSpecies)

#Beta dis calculated between sites within a year
#source('saved/results/partnerVar_Site.Rdata')

#much more variable in 2012 than later years?
beta.dist %>% 
  filter(Year == 2022) %>% 
  ggplot(
    aes(x = dist)) +
  geom_histogram() +
  facet_wrap( ~ Site)

beta.dist %>% 
  ggplot(
    aes(x = Site, y = dist)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Year)

print(paste("Plant species", length(unique(beta.dist$GenusSpecies)))) #48 species
unique(beta.dist$GenusSpecies)

