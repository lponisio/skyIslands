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

load('data/splevel_network_metrics/YearSR_PlantPollinator_Bees.Rdata')


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

pol.pca.scores[1]


autoplot(plant.pca.scores$'2018'$pca.loadings, loadings=TRUE,
         loadings.colour = 'blue')


save(plant.pca.scores,  file="analysis/role/saved/results/pol_pcaVar.Rdata")

##################################################################
## --- Merge role variability, network metrics, and climate --- ##
##################################################################
load('analysis/role/saved/results/beeSyrphid_pcaVar.Rdata')
load('data/spec_traits.Rdata')
getwd()

#combine PCA results from all years into a single dataframe
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
  mutate(beeEmergence_start = min(Doy),
         beeEmergence_end = max(Doy),
         flight_period = beeEmergence_end - beeEmergence_start,
         .groups = "drop") %>% 
  group_by(GenusSpecies, Site, Year) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), 
            .groups = "drop") %>% 
  mutate(GenusSpeciesSiteYear = paste0(GenusSpecies, Site, Year))

## --- merge species level network metrics, partner variability, trait, and climate data --- ##
pol.pcas.network <- spec %>%
  right_join(pol.pcas, by = "GenusSpeciesSiteYear") # keeps all rows from pol.beta
colnames(pol.pcas.network)

write.csv(pol.pcas.network, file = 'analysis/role/saved/traits/pol_pcas.csv')

## --- merge climate data --- ##
climate <- read.csv('analysis/role/saved/traits/climateVariability.csv')

pol.pcas.climate <- pol.pcas.network %>%
  left_join(climate, by = c("Site", "Year"))

write.csv(pol.pcas.climate, file = 'saved/traits/pol_pcas_climate.csv')








