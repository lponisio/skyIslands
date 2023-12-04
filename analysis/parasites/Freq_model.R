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

## Freq function by Rebecca
run_plot_freq_model_diagnostics <- function(this_formula, #brms model formula
                                            this_data, #data frame, subsetted to correct weights!
                                            num_chains=1,
                                            num_iter=10000,
                                            this_family #model family
){
  #run model
  this_model_output <- brms::brm(this_formula,
                                 data = this_data,
                                 chains = num_chains,
                                 iter = num_iter, family=this_family)
  this_model_output
  # return a list of single plots
  diagnostic.plots <- plot(check_model(this_model_output, panel = TRUE))
  diagnostic.plots
}
## Frequentist formulas
## **********************************************************
## Model 1.1: formula for lat effects on floral community
## **********************************************************
## flower diversity
freq.flower.div <- formula(MeanFloralDiversity ~
                                Lat +
                                Year +
                                SRDoy + I(SRDoy^2) +
                                (1|Site)
)
## flower abund
freq.flower.abund <- formula(MeanFloralAbundance ~
                                  Year+ Lat + 
                                  SRDoy + I(SRDoy^2) +
                                  (1|Site)
)

## **********************************************************
## Model 1.2: formula for forest effects on bee community
## **********************************************************
## bee diversity
freq.bee.div <- formula(Net_BeeDiversity ~
                             MeanFloralDiversity +
                             Lat + Year +
                             SRDoy + I(SRDoy^2) +
                             (1|Site)
)
## bombus abund
freq.bombus.abund <- formula(Net_BombusAbundance ~
                                  MeanFloralAbundance + Year +
                                  SRDoy + I(SRDoy^2) +
                                  Lat + 
                                  (1|Site)
)
## HB abund
freq.HB.abund <- formula(Net_HBAbundance ~
                              MeanFloralAbundance +  Year +
                              SRDoy + I(SRDoy^2) +
                              Lat +
                              (1|Site)
)
## bee abund
freq.bee.abund <- formula(Net_NonBombusHBAbundance ~
                               MeanFloralAbundance +  Year +
                               SRDoy + I(SRDoy^2) +
                               Lat +
                               (1|Site)
)

## **********************************************************
## convert formulas to brms forma
## **********************************************************
bf.fabund <- bf(formula.flower.abund)
bf.fdiv <- bf(formula.flower.div)
bf.babund <- bf(formula.bee.abund, family="negbinomial")
bf.bombusabund <- bf(formula.bombus.abund, family="negbinomial")
bf.HBabund <- bf(formula.HB.abund, family="negbinomial")
bf.bdiv <- bf(formula.bee.div)

form.community <- bf.fabund + bf.fdiv +
  bf.babund +
  bf.bombusabund + bf.HBabund +
  bf.bdiv  +
  set_rescor(FALSE)

freq.community.model <- run_plot_freq_model_diagnostics(form.community,
                                                         spec.net[spec.net$Weights==1,], this_family = 'gaussian')






## Write the parasite formula using the beta binomial

beta_bi_formula <- formula(CrithidiaPresence|trials(1) ~ Net_NonBombusHBAbundance + 
                             Net_HBAbundance + Net_BombusAbundance + 
                             Net_BeeDiversity + rare.degree + MeanITD + 
                             (1|Site) + (1|GenusSpecies)) 
