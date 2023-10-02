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
vars <- c("MeanFloralAbundance",
          "MeanFloralDiversity",
          "Net_BeeDiversity",
          "Lat", "SRDoy",
          "MeanITD",
          ## "r.degree", ## across all networks
          "rare.degree"  ## site-year level
)

variables.to.log <- c("MeanFloralAbundance",
                      "Net_NonBombusHBAbundance",
                      "Net_HBAbundance",
                      "Net_BombusAbundance",
                      "Lat")
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
hist(spec.net$Net_BeeAbundance)
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
##funciton to run frequentist version of brms models and plot diagnostics
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

#Floral Diversity brms formula
formula.flower.div <- formula(MeanFloralDiversity | weights(Weights) ~
                                Lat + SRDoy + I(SRDoy^2) + Year +
                                (1|Site)
)
#Floral Diversity freq formula
freq.formula.flower.div <- formula(MeanFloralDiversity ~
                                    Lat + SRDoy + I(SRDoy^2) + Year +
                                    (1|Site)
)
#for this_data, use spec.net[spec.net$Weights==1,] to incorporate weights into frequentist models
freq.flower.div.model <- run_plot_freq_model_diagnostics(freq.formula.flower.div,
                                                         spec.net[spec.net$Weights==1,], this_family = 'gaussian')
#dir.create("figures", showWarnings = FALSE)
ggsave(freq.flower.div.model, file="figures/FloralDiversityModelDiagnostics.pdf",
       height=8, width=11)

## flower abund
formula.flower.abund <- formula(MeanFloralAbundance  ~
                                  SRDoy + I(SRDoy^2) + Year +
                                  (1|Site)
)

freq.flower.abun.model <- run_plot_freq_model_diagnostics(formula.flower.abund,
                                                         spec.net[spec.net$Weights==1,], this_family = 'gaussian')
ggsave(freq.flower.abun.model, file="figures/FloralAbunModelDiagnostics.pdf",
       height=8, width=11)

## bee diversity
formula.bee.div <- formula(Net_BeeDiversity ~
                             MeanFloralDiversity +
                             Lat +  SRDoy + Year +
                             (1|Site)
)
freq.bee.div.model <- run_plot_freq_model_diagnostics(formula.bee.div,
                                                          spec.net[spec.net$Weights==1,], this_family = 'gaussian')
ggsave(freq.bee.div.model, file="figures/BeeDivModelDiagnostics.pdf",
       height=8, width=11)

## bombus abund
formula.bombus.abund <- formula(Net_BombusAbundance ~
                                  MeanFloralAbundance +
                                  MeanFloralDiversity+
                                  SRDoy + I(SRDoy^2) +
                                  Lat + Year+
                                  (1|Site)
)
freq.bombus.abun.model <- run_plot_freq_model_diagnostics(formula.bombus.abund,
                                                          spec.net[spec.net$Weights==1,], this_family = 'gaussian')
ggsave(freq.bombus.abun.model, file="figures/BombusAbunModelDiagnostics.pdf",
       height=8, width=11)

## HB abund
formula.HB.abund <- formula(Net_HBAbundance ~
                              MeanFloralAbundance +
                              MeanFloralDiversity+
                              SRDoy + I(SRDoy^2) +
                              Lat + Year+
                              (1|Site)
)
freq.HB.abun.model <- run_plot_freq_model_diagnostics(formula.HB.abund,
                                                          spec.net[spec.net$Weights==1,], this_family = 'gaussian')
ggsave(freq.HB.abun.model, file="figures/HBAbunModelDiagnostics.pdf",
       height=8, width=11)

## bee abund
formula.bee.abund <- formula(Net_NonBombusHBAbundance ~
                               MeanFloralAbundance +
                               MeanFloralDiversity+
                               SRDoy + I(SRDoy^2) +
                               Lat + Year+
                               (1|Site)
)
freq.bee.abun.model <- run_plot_freq_model_diagnostics(formula.bee.abund,
                                                      spec.net[spec.net$Weights==1,], this_family = 'gaussian')
ggsave(freq.bee.abun.model, file="figures/NonBombusHBAbunModelDiagnostics.pdf",
       height=8, width=11)

## parasite
formula.parasite <- formula(parasites ~
                                MeanFloralAbundance+
                              MeanFloralDiversity+
                              Net_BeeDiversity+
                              Lat+ SRDoy+
                              MeanITD +
                              r.degree +
                              rare.degree 
)
freq.parasite.model <- run_plot_freq_model_diagnostics(formula.parasite,
                                                       spec.net[spec.net$Weights==1,], this_family = 'gaussian')
ggsave(freq.parasite.model, file="figures/ParasiteModelDiagnostics.pdf",
       height=8, width=11)