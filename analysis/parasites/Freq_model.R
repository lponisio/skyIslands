rm(list=ls())
setwd("C:/Users/na_ma/Dropbox (University of Oregon)/skyIslands")
library(loo)
library(brms)
library(bayesplot)
library(rstanarm)
library(lme4)
library(performance)

setwd("analysis/parasites")
source("src/misc.R")
source("src/writeResultsTable.R")
source("src/makeMultiLevelData.R")
source("src/runParasiteModels.R")
source("src/standardize_weights.R")
## all of the variables that are explanatory variables and thus need
## to be centered
vars_yearsr <- c("MeanFloralAbundance",
                 "MeanFloralDiversity",
                 "Net_BeeDiversity",
                 "Lat", "SRDoy"  
)
vars_sp <- c("MeanITD",
             "rare.degree")

variables.to.log <- "rare.degree"

## uses only net specimens, and drops syrphids
source("src/init.R")

