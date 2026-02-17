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


## calculate emergence and flight period ##
load('saved/results/pol_partnerVar_Site.Rdata')
load('../../data/spec_traits.Rdata')

spec <- spec.net %>% 
  group_by(GenusSpecies, Year) %>% 
  mutate(beeEmergence_start = min(Doy),
         beeEmergence_end = max(Doy),
         flight_period = beeEmergence_end - beeEmergence_start,
         .groups = "drop") %>% 
  mutate(GenusSpeciesSiteYear = paste0(GenusSpecies, Site, Year))

pol.partner <- beta.dist %>%
  mutate(GenusSpeciesSiteYear = paste0(GenusSpecies, Site, Year))

## merge species level network metrics, role variability, trait data ##
pol.beta.traits <- spec %>%
  right_join(pol.partner, by = "GenusSpeciesSiteYear") # keeps all rows from pol.beta

write.csv(pol.beta.traits, file = 'saved/traits/pol_partner_traits.csv')


## create vector with columns and filter unwanted columns ##
colnames(pol.beta.traits)

names <- c('')