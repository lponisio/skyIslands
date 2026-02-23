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
load('../../data/spec_traits_year.Rdata')

spec <- spec.net %>% 
  group_by(GenusSpecies, Year) %>% 
  mutate(PolEmergenceStart = min(Doy),
         PolEmergenceEnd = max(Doy),
         FlightPeriod = PolEmergenceEnd - PolEmergenceStart) %>%
  ungroup() %>% 
  mutate(GenusSpeciesSiteYear = paste0(GenusSpecies, Site, Year))

pol.partner <- beta.dist %>%
  mutate(GenusSpeciesSiteYear = paste0(GenusSpecies, Site, Year))

## merge species level network metrics, role variability, trait data ##
pol.beta.traits <- spec %>%
  right_join(pol.partner %>% select(-SampleRound, -Site, -Year, -GenusSpecies), 
             by = "GenusSpeciesSiteYear")

## create vector with columns names and filter unwanted columns ##
names <- c('SiteYearSR', 'SampleRound', 'SiteYear', 'Site', 'Year', 'Date', 'Lat', 'Long', 'Elev', 'Area', 'GenusSpecies',
           'Sex', 'Order', 'Family', 'Genus', 'SubGenus', 'Species', 'SubSpecies',
           'PlantGenus', 'PlantSpecies', 'PlantVar', 'PlantSubSpecies', 'State', 'County', 'Meadow', 'Forest',
           'MtRange', 'PlantGenusSpecies', 'Int', 'Doy', 'PlantFamily', 'degree',
           'normalised.degree', 'species.strength', 'interaction.push.pull', 'nestedrank',
           'PDI', 'resource.range', 'species.specificity.index', 'PSI', 'node.specialisation.index.NSI',
           'betweenness', 'weighted.betweenness', 'closeness', 'weighted.closeness', 'Fisher.alpha',
           'partner.diversity', 'effective.partners', 'proportional.generality', 'proportional.similarity', 
           'd', 'tot.int', 'niche.overlap', 'rare.degree', 'PollAbundance', 'BeeAbundance', 'SyrphidAbundance', 
           'HBAbundance', 'BombusAbundance', 'NonBombusHBAbundance', 'PollRichness', 'BeeRichness',
           'SyrphidRichness', 'BombusRichness', 'PollDiversity', 'BeeDiversity', 'SyrphidDiversity',
           'BombusDiversity', 'VisitedFloralRichness', 'VisitedFloralDiversity', 'MeanFloralRichness',
           'MeanFloralDiversity', 'MeanFloweringPlantDiversity', 'MeanFloweringPlantAbundance', 'Abundance', 'uniqSite',
           'originalitySite', 'originalitySiteYearSR', 'BeeFDisSiteYearSR', 'BeeFEveSiteYearSR', 'uniqSiteYear',
           'originalitySiteYear', 'BeeFDisSiteYear', 'BeeFEveSiteYear', 'PolEmergenceStart', 'PolEmergenceEnd',
           'FlightPeriod', 'GenusSpeciesSiteYear', 'dist', 'Type')

pol.beta.traits <- pol.beta.traits[, names]

write.csv(pol.beta.traits, file = 'saved/traits/Bee_partner.csv')



## Make dataset for bee and syrphids ##
load('saved/results/pol_partnerVar_Site.Rdata')
load('../../data/spec_net.Rdata')

## calculate flight period ##
spec <- spec.net %>% 
  group_by(GenusSpecies, Year) %>% 
  mutate(PolEmergenceStart = min(Doy),
         PolEmergenceEnd = max(Doy),
         FlightPeriod = PolEmergenceEnd - PolEmergenceStart) %>%
  ungroup() %>% 
  mutate(GenusSpeciesSiteYear = paste0(GenusSpecies, Site, Year))

pol.partner <- beta.dist %>%
  mutate(GenusSpeciesSiteYear = paste0(GenusSpecies, Site, Year))

## merge species level network metrics, role variability, trait data ##
pol.beta.traits <- spec %>%
  right_join(pol.partner %>% select(-SampleRound, -Site, -Year, -GenusSpecies), 
             by = "GenusSpeciesSiteYear")

bee.syrphid.familes <- c("Halictidae", "Syrphidae", "Apidae", "Bombyliidae", "Andrenidae", "Megachilidae", "Colletidae")

pol.beta.traits <- pol.beta.traits %>% 
  filter(Family %in% bee.syrphid.familes)

unique(pol.beta.traits$Family)

## create vector with columns names and filter unwanted columns ##
names <- c('SampleRound', 'Site', 'Year', 'Date', 'Lat', 'Long', 'Elev', 'Area', 'GenusSpecies',
           'Sex', 'Order', 'Family', 'Genus', 'SubGenus', 'Species', 'SubSpecies',
           'PlantGenus', 'PlantSpecies', 'PlantVar', 'PlantSubSpecies', 'State', 'County', 'Meadow', 'Forest',
           'MtRange', 'PlantGenusSpecies', 'Int', 'Doy', 'PlantFamily', 'degree',
           'normalised.degree', 'species.strength', 'interaction.push.pull', 'nestedrank',
           'PDI', 'resource.range', 'species.specificity.index', 'PSI', 'node.specialisation.index.NSI',
           'betweenness', 'weighted.betweenness', 'closeness', 'weighted.closeness', 'Fisher.alpha',
           'partner.diversity', 'effective.partners', 'proportional.generality', 'proportional.similarity', 
           'd', 'tot.int', 'niche.overlap', 'rare.degree', 'PollAbundance', 'BeeAbundance', 'SyrphidAbundance', 
           'HBAbundance', 'BombusAbundance', 'NonBombusHBAbundance', 'PollRichness', 'BeeRichness',
           'SyrphidRichness', 'BombusRichness', 'PollDiversity', 'BeeDiversity', 'SyrphidDiversity',
           'BombusDiversity', 'VisitedFloralRichness', 'VisitedFloralDiversity', 'MeanFloralRichness',
           'MeanFloralDiversity', 'MeanFloweringPlantDiversity', 'MeanFloweringPlantAbundance', 'Abundance', 'PolEmergenceStart', 'PolEmergenceEnd',
           'FlightPeriod', 'GenusSpeciesSiteYear', 'dist', 'Type')

pol.beta.traits <- pol.beta.traits[, names]

write.csv(pol.beta.traits, file = 'saved/traits/BeeSyrphid_partner.csv')
