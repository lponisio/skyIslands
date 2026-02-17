## This script calculates the pollinator role/network niche
## variability between sites within a year.
rm(list=ls())
library(ggfortify)
library(bipartite)
library(fossil)
setwd("C:/")
source("lab_paths.R")
local.path

dir.bombus <- file.path(local.path, "skyIslands")
setwd(dir.bombus)

this.script <- "role"
source('analysis/role/src/initialize.R')
type <- "all"


## vector of pca loadings of interest
loadings <- c(1)
metrics <- c("rare.degree",
             "weighted.betweenness",
             "weighted.closeness",
             "niche.overlap",
             "species.strength",
             "d")

## Metrics used in the PCA
var.method <- cv
ave.method <- mean


## PCA 
pol.pca.scores <- calcPcaMeanVar(species.roles=sp.lev, 
                                 var.method=var.method,
                                 ave.method=ave.method,
                                 metrics= metrics,
                                 loadings=loadings,
                                 agg.col = "Year")


autoplot(plant.pca.scores$'2018'$pca.loadings, loadings=TRUE,
         loadings.colour = 'blue')


save(pol.pca.scores,  file="analysis/role/saved/results/beeSyrphid_pcaVar.Rdata")

##################################################################
## --- Merge role variability, network metrics, and climate --- ##
##################################################################
load('analysis/role/saved/results/beeSyrphid_pcaVar.Rdata')
load('data/spec_traits.Rdata')

## combine PCA results from all years into a single dataframe ##
pol.pcas <- do.call( #combines results and applies function
  rbind, #row bind list of dataframes returned by lapply()
  lapply(names(pol.pca.scores), # Loop over each list element name (year "2012")
         function(yr) { 
           df <- pol.pca.scores[[yr]]$pcas #extract the PCA scores dataframe for this year
           df$Year <- yr #add a column recording which list element (year)
           return(df)}))

pol.pcas <- pol.pcas %>%
  mutate(GenusSpeciesSiteYear = paste0(GenusSpecies, Site, Year))

spec <- spec.net %>% 
  group_by(GenusSpecies, Year) %>% 
  mutate(PolEmergenceStart = min(Doy),
         PolEmergenceEnd = max(Doy),
         FlightPeriod = PolEmergenceEnd - PolEmergenceStart,
         .groups = "drop") %>% 
  # group_by(GenusSpecies, Site, Year) %>%
  # summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), 
  #           .groups = "drop") %>% 
  mutate(GenusSpeciesSiteYear = paste0(GenusSpecies, Site, Year))

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
           'FlightPeriod', 'GenusSpeciesSiteYear')

spec <- spec[, names]

## merge species level network metrics, partner variability, trait, and climate data ##
pol.pcas.network <- spec %>%
  right_join(pol.pcas %>%  select(-Site, -Year, -GenusSpecies),
             by = "GenusSpeciesSiteYear") 

colnames(pol.pcas.network)

write.csv(pol.pcas.network, file = 'analysis/role/saved/traits/BeeSyrphid_pcas.csv')

## Make dataset for only bee familes ##
bee.families <- c("Andrenidae", "Apidae", "Colletidae", "Halictidae",
                  "Megachilidae")

bee.pcas.network <- pol.pcas.network %>% 
  filter(Family %in% bee.families)

write.csv(bee.pcas.network, file = 'analysis/role/saved/traits/Bee_pcas.csv')








