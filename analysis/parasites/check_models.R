rm(list=ls())
setwd("C:/Users/na_ma/Dropbox (University of Oregon)/skyIslands")
library(loo)
library(brms)
library(bayesplot)
library(rstanarm)
library(lme4)
library(performance)
library(glmmTMB)

setwd("analysis/parasites")
source("src/misc.R")
source("src/writeResultsTable.R")
source("src/makeMultiLevelData.R")
source("src/runParasiteModels.R")
source("src/standardize_weights.R")
source("src/runPlotFreqModelDiagnostics.R")
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
variables.to.log.1 <- c("Net_HBAbundance", "Net_BombusAbundance", 
                        "Net_NonBombusHBAbundance")

## uses only net specimens, and drops syrphids
source("src/init.R")
## Testing the fit of the parasite models. 


## Models
## **********************************************************
## Model 1.1: formula for forest effects on floral community
## **********************************************************
## flower diversity
formula.flower.div <- formula(MeanFloralDiversity | weights(Weights) ~
                                Lat + SRDoy + I(SRDoy^2) + Year +
                                (1|Site)
)

#checking posterior
flowercheck1 <- brms::brm(formula.flower.div,
                          data = spec.net, chains = 1,
                          iter = 1000)
mcmc_areas(as.matrix(flowercheck1),
           prob_outer = .999)
##LOO validation - checking pareto K plot
flowerloop <- loo(flowercheck1)
plot(flowerloop)

## flower abund
formula.flower.abund <- formula(log(MeanFloralAbundance) | weights(Weights) ~
                                  SRDoy + I(SRDoy^2) + Year +
                                  (1|Site)
)
#checking posterior
flowercheck2 <- brms::brm(formula.flower.div,
                          data = spec.all, chains = 1,
                          iter = 1000)
mcmc_areas(as.matrix(flowercheck2),
           prob_outer = .999)
##LOO validation - checking pareto K plot
flowerloop <- loo(flowercheck2)
plot(flowerloop)

################################################################################
## Histogram of raw data
hist(spec.net$MeanFloralDiversity)
hist(spec.net$MeanFloralAbundance)
hist(log(spec.net$Net_BombusAbundance)+ 1)
hist(log(spec.net$Net_HBAbundance) + 1)
hist(spec.net$Net_BeeDiversity)

spec.net$MeanFloralDiversity<- log(spec.net$MeanFloralDiversity)
hist(spec.net$MeanFloralDiversity)

spec.net$MeanFloralAbundance<- log(spec.net$MeanFloralAbundance)
hist(spec.net$MeanFloralAbundance)




formula.flower.div <- formula(MeanFloralDiversity | weights(Weights) ~
                                Lat + SRDoy + I(SRDoy^2) + Year +
                                (1|Site)
)

#checking posterior
flowercheck1 <- brms::brm(formula.flower.div,
                          data = spec.net, chains = 1,
                          iter = 1000)
## Get residuals
spec.net %>%
  add_residual_draws(flowercheck1) %>%
  ggplot(aes(x = .row, y = .residual)) +
  stat_pointinterval()


spec.net %>%
  add_residual_draws(flowercheck1) %>%
  median_qi() %>%
  ggplot(aes(sample = .residual)) +
  geom_qq() +
  geom_qq_line()  

###############################################################################
## Using the performance package to get the R2
## Function by Rebecca


#Floral Diversity brms formula
formula.flower.div <- formula(MeanFloralDiversity | weights(Weights) ~
                                Lat +
                                  SRDoy*Year + I(SRDoy^2)*Year +
                                  SRDoy*Lat + I(SRDoy^2)*Lat +
                                  (1|Site)
)
#Floral Diversity freq formula
formula.flower.div <- formula(MeanFloralDiversity ~
                                Lat +
                                Year +
                                SRDoy + I(SRDoy^2) +
                                (1|Site)
)
#for this_data, use spec.net[spec.net$Weights==1,] to incorporate weights into frequentist models
freq.flower.div.model <- run_plot_freq_model_diagnostics(freq.formula.flower.div,
                                                         spec.net[spec.net$Weights==1,], this_family = 'gaussian')
#dir.create("figures", showWarnings = FALSE)
ggsave(freq.flower.div.model, file="figures/FloralDiversityModelDiagnostics.pdf",
       height=8, width=11)

## flower abund
formula.flower.abund <- formula(MeanFloralAbundance ~
                                                          Year+ Lat + 
                                                          SRDoy + I(SRDoy^2) +
                                                          (1|Site)
)

freq.flower.abun.model <- run_plot_freq_model_diagnostics(formula.flower.abund,
                                                         spec.net[spec.net$Weights==1,], this_family = 'gaussian')
ggsave(freq.flower.abun.model, file="figures/FloralAbunModelDiagnostics.pdf",
       height=8, width=11)

## bee diversity
formula.bee.div <- formula(Net_BeeDiversity ~
                             MeanFloralDiversity +
                             Lat + Year +
                             SRDoy + I(SRDoy^2) +
                             (1|Site)
)
freq.bee.div.model <- run_plot_freq_model_diagnostics(formula.bee.div,
                                                          spec.net[spec.net$Weights==1,], this_family = 'gaussian')
ggsave(freq.bee.div.model, file="figures/BeeDivModelDiagnostics.pdf",
       height=8, width=11)

## bombus abund
formula.bombus.abund <- formula(Net_BombusAbundance ~
                                  MeanFloralAbundance + 
                                  #Year +
                                  SRDoy + #I(SRDoy^2) +
                                  Lat + 
                                  (1|Site)
)
freq.bombus.abun.model <- run_plot_freq_model_diagnostics(formula.bombus.abund,
                                                          spec.net[spec.net$Weights==1,], this_family = 'gaussian')
ggsave(freq.bombus.abun.model, file="figures/BombusAbunModelDiagnostics_gaussian.pdf",
       height=8, width=11)

## HB abund
formula.HB.abund <- formula(Net_HBAbundance ~
                              MeanFloralAbundance +  
                              SRDoy + #I(SRDoy^2) +
                              Lat +
                              (1|Site)
)
freq.HB.abun.model <- run_plot_freq_model_diagnostics(formula.HB.abund,
                                                          spec.net[spec.net$Weights==1,], this_family = "gaussian")
ggsave(freq.HB.abun.model, file="figures/HBAbunModelDiagnostics_gaussian.pdf",
       height=8, width=11)

## bee abund
formula.bee.abund <- formula(Net_NonBombusHBAbundance ~
                               MeanFloralAbundance +
                               SRDoy  +
                               Lat +
                               (1|Site)
)
freq.bee.abun.model <- run_plot_freq_model_diagnostics(formula.bee.abund,
                                                      spec.net[spec.net$Weights==1,], this_family = 'gaussian')
ggsave(freq.bee.abun.model, file="figures/NonBombusHBAbunModelDiagnostics_gaussian.pdf",
       height=8, width=11)

## parasite
beta_bi_formula <- formula(CrithidiaPresence ~ Net_NonBombusHBAbundance + 
                             Net_HBAbundance + Net_BombusAbundance + 
                             Net_BeeDiversity + rare.degree + MeanITD + 
                             (1|Site) + (1|GenusSpecies)) 

freq.parasite.model <- run_plot_freq_model_diagnostics(beta_bi_formula,
                                                       spec.net[spec.net$WeightsPar==1,], 
                                                       this_family = 'bernoulli')

ggsave(freq.parasite.model, file="figures/ParasiteModelDiagnostics.pdf",
       height=8, width=11)


## Load the model data
load(file="saved/communityFit.Rdata")
## Posterior predictive check for the community model with students t distribution
pp_check(fit.community, Net_BombusAbundance)
