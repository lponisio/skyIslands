rm(list=ls())
## setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')
setwd("C:/Users/na_ma/Dropbox (University of Oregon)/skyIslands")
##setwd('~/Dropbox (University of Oregon)/skyislands')
library(brms)
library(tidybayes)
library(stringr)
setwd("analysis/parasites")
load(file="saved/spec_weights.Rdata")
source("src/misc.R")

load("saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat_all_bees.Rdata")
bombus.all <- fit.parasite
load("saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat_social_species.Rdata")
bombus.ss <- fit.parasite
load("saved/parasiteFit_apis_CrithidiaPresenceApicystisSpp_lat_apis_abundance.Rdata")
apis.ha <- fit.parasite

## Generate newdata draws

bombus.all.cond.effects <- conditional_effects(bombus.all)
bombus.ss.cond.effects <- conditional_effects(bombus.ss)
apis.cond.effects <- conditional_effects(apis.ha)


## create a function that takes the percentage difference between the lowest and
## the highest point

CalcDiff <- function(cond_effect, effect = "char"){
  # This function calculates the percentage difference between the lowest and 
  # highest points
  #cond_effec = list of models with the "fake" data created from the 
  #posterior of the brms model
  #effect = character specifying the parameter in list cond_effect
  parameter_obs <- cond_effect[[effect]] 
  V1 <- min(parameter_obs$estimate__, na.rm = TRUE)
  V2 <- max(parameter_obs$estimate__, na.rm = TRUE)
  diff <- abs(V1 - V2)
  return(list(effect = effect, diff = diff))
}

## Vector with all the names of the parameters in the parasite models
par_predictors <- c("CrithidiaPresence.CrithidiaPresence_MeanFloralDiversity",
                    "ApicystisSpp.ApicystisSpp_MeanFloralDiversity",
                    "CrithidiaPresence.CrithidiaPresence_Net_BeeDiversity",
                    "ApicystisSpp.ApicystisSpp_Net_BeeDiversity",
                    "CrithidiaPresence.CrithidiaPresence_Net_BeeAbundance",
                    "ApicystisSpp.ApicystisSpp_Net_BeeAbundance",
                    "CrithidiaPresence.CrithidiaPresence_rare.degree",
                    "ApicystisSpp.ApicystisSpp_rare.degree",
                    "CrithidiaPresence.CrithidiaPresence_Lat",
                    "ApicystisSpp.ApicystisSpp_Lat"
)

## Empty list to put the values calculated
bombus_percentages <- vector("list", length = length(par_predictors))
## For loop that calculates the percentage difference through all the parasite parameters
for (i in seq_along(par_predictors)) {
    bombus_percentages[[i]] <- 
      CalcDiff(bombus.all.cond.effects, effect = par_predictors[i])
}

## Repeat: For Apicystis model
## Vector with all the names of the parameters in the parasite models
par_predictors <- c("CrithidiaPresence.CrithidiaPresence_MeanFloralDiversity",
                    "ApicystisSpp.ApicystisSpp_MeanFloralDiversity",
                    "CrithidiaPresence.CrithidiaPresence_Net_BeeDiversity",
                    "ApicystisSpp.ApicystisSpp_Net_BeeDiversity",
                    "CrithidiaPresence.CrithidiaPresence_Net_BombusAbundance",
                    "CrithidiaPresence.CrithidiaPresence_Net_HBAbundance",
                    "ApicystisSpp.ApicystisSpp_Net_BombusAbundance",
                    "ApicystisSpp.ApicystisSpp_Net_HBAbundance",
                    "CrithidiaPresence.CrithidiaPresence_rare.degree",
                    "ApicystisSpp.ApicystisSpp_rare.degree",
                    "CrithidiaPresence.CrithidiaPresence_Lat",
                    "ApicystisSpp.ApicystisSpp_Lat"
)

## Empty list to put the values calculated
bapicystis_percentages <- vector("list", length = length(par_predictors))
## For loop that calculates the percentage difference through all the parasite parameters
for (i in seq_along(par_predictors)) {
  bapicystis_percentages[[i]] <- 
    CalcDiff(bombus.ss.cond.effects, effect = par_predictors[i])
}

## Repeat: Calculate for apis models
apis_par_predictors <- c("CrithidiaPresence.CrithidiaPresence_Net_BeeDiversity",
                    "ApicystisSpp.ApicystisSpp_Net_BeeDiversity",
                    "CrithidiaPresence.CrithidiaPresence_Net_HBAbundance",
                    "ApicystisSpp.ApicystisSpp_Net_HBAbundance",
                    "CrithidiaPresence.CrithidiaPresence_rare.degree",
                    "ApicystisSpp.ApicystisSpp_rare.degree",
                    "CrithidiaPresence.CrithidiaPresence_Lat",
                    "ApicystisSpp.ApicystisSpp_Lat"
)
apis_percentages <- vector("list", length = length(apis_par_predictors))

for (i in seq_along(apis_par_predictors)) {
  apis_percentages[[i]] <- 
    CalcDiff(apis.cond.effects, effect = apis_par_predictors[i])
}
