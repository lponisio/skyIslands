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

#combine PCA results from all years into a single dataframe
all.pcas <- do.call( #combines results and applies function
  rbind, #row bind list of dataframes returned by lapply()
  lapply(names(pol.pca.scores), # Loop over each list element name (year "2012")
         function(yr) { 
    df <- pol.pca.scores[[yr]]$pcas #extract the PCA scores dataframe for this year
    df$Year <- yr #add a column recording which list element (year)
    return(df)})) #returns dataframe

## *********************************************************
## Create summary figures
## *********************************************************

library(ggplot2)
library(tidyverse)

## Plants

all.pcas %>%
  filter(!is.na(var.pca1)) %>% 
  ggplot(aes(x = var.pca1)) +
    geom_histogram()


outliers <- all.pcas %>% 
  filter(var.pca1 > 200 | var.pca1 < -20)

var.pcas <- all.pcas %>% 
  filter(!(var.pca1 > 200 | var.pca1 < -20))

#Variance
var.pcas %>% 
  ggplot(aes(x = var.pca1)) +
  geom_histogram() +
  facet_wrap(~Site)

var.pcas %>% 
  ggplot(aes(x = var.pca1)) +
  geom_histogram() +
  facet_wrap(~Year)

var.pcas %>% 
  arrange(var.pca1) %>% 
  ggplot(
    aes(y = GenusSpecies, x = var.pca1)) +
  geom_point() +
  labs(x = "Partner variability", y = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


print(paste("Plant species", length(unique(all.pcas$GenusSpecies)))) #100 plant species
unique(all.pcas$GenusSpecies)

## Pollinators

var.pcas %>%
  filter(!is.na(var.pca1)) %>% 
  ggplot(aes(x = var.pca1)) +
  geom_histogram()

var.pcas <- all.pcas %>% 
  filter(!(var.pca1 < -100))

var.pcas %>% 
  ggplot(
    aes(y = reorder(GenusSpecies, var.pca1), x = var.pca1)) +
  geom_point() +
  labs(x = "Partner variability", y = "Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
