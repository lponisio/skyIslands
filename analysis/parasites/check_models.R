rm(list=ls())
setwd("C:/Users/na_ma/Dropbox (University of Oregon)/skyIslands")
library(loo)
library(brms)
library(bayesplot)
library(rstanarm)

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
library(performance)


formula.flower.div <- formula(MeanFloralDiversity | weights(Weights) ~
                                Lat + SRDoy + I(SRDoy^2) + Year +
                                (1|Site)
)
library(lme4)
model <- lmer(MeanFloralDiversity  ~
                Lat + SRDoy + I(SRDoy^2) + Year +
                (1|Site), data = spec.net)
r2(model)

model_performance(model)
check_model(model)




