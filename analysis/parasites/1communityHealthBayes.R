setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')

## setwd('~/Dropbox (University of Oregon)/skyislands')

setwd("analysis/parasites")
rm(list=ls())

source("src/init.R")
source("src/writeResultsTable.R")
source("src/makeMultiLevelData.R")
source("src/misc.R")

ncores <- 10

## **********************************************************
## formula for site effects on the bee community
## **********************************************************

## order data by stand for ease of viewing
spec <- spec[order(spec$Site),]


## all of the variables that are explanatory variables and thus need
## to be centered
vars <- c("FloralAbundance",
          "FloralDiversity",
          "PollAbundance",
          "PollDiversity",
          "Elev",
          "Lat")

##  center all of the x variables across the datasets
spec[, vars] <- apply(spec[, vars], 2, standardize)

## create a dumby varaible "Weight" to deal with the data sets being at
## different levels to get around the issue of having to pass in one
## data set into brms

## will need to modify when we have multiple years
spec <- makeDataMultiLevel(spec, "Site", "Year")

## create a dumby varaible "WeightPar" for the parasite data. The
## original intention was to keep stan from dropping data for
## site-level models, but weight is 0 for parasite models.

spec$WeightsPar <- 1
spec$WeightsPar[is.na(spec$ParasitePresence) | spec$Apidae != 1] <- 0

## stan drops all NA data, so can set ParasitePresence to 0 with WeightsPar
## to keep it in the models, this is commented out because we don't
## want to use the entire dataset in this publication

spec$ParasitePresence[is.na(spec$ParasitePresence | spec$Apidae != 1)] <- 0

## define all the formulas for the different parts of the models
## **********************************************************
## Model 1.1: formula for forest effects on floral community
## **********************************************************
## flower diversity
formula.flower.div <- formula(FloralDiversity | weights(Weights) ~
                                  Elev + Lat
                              )
## flower abund
formula.flower.abund <- formula(FloralAbundance | weights(Weights) ~
                                    Lat + Elev
                                )

## **********************************************************
## Model 1.2: formula for forest effects on bee community
## **********************************************************
## bee diversity
formula.bee.div <- formula(PollDiversity | weights(Weights)~
                               FloralAbundance +
                               FloralDiversity +
                               Lat)



## bee abund
formula.bee.abund <- formula(PollAbundance | weights(Weights)~
                                 FloralAbundance +
                                 FloralDiversity +
                                 Lat)

## **********************************************************
## Model 1.3: formula for bee community effects on parasitism
## **********************************************************
formula.parasite <- formula(ParasitePresence | weights(WeightsPar) ~
                                PollAbundance +
                                FloralDiversity +
                                PollDiversity +
                                FloralAbundance +
                                (1|Site)
                            )

## **********************************************************
## Community models
## **********************************************************
## convert to brms format
bf.fabund <- bf(formula.flower.abund)
bf.fdiv <- bf(formula.flower.div)
bf.babund <- bf(formula.bee.abund)
bf.bdiv <- bf(formula.bee.div)

## **********************************************************
## Model 1 community effects on bee parasitism
## **********************************************************
bf.par <- bf(formula.parasite, family="bernoulli")

## full model
bform <- bf.fabund + bf.fdiv + bf.babund + bf.bdiv + bf.par +
    set_rescor(FALSE)

## run model
fit <- brm(bform, spec,
           cores=ncores,
           iter = 10^4,
           chains = 2,
           thin=1,
           inits=0,
           control = list(adapt_delta = 0.99))

write.ms.table(fit, "parasitism")

save(fit, spec,
     file="saved/parasiteFitMod.Rdata")

## dignostic figures
mcmc_trace(fit)
ggsave("figures/diagnostics/parasite.pdf",
       height=11, width=8.5)


## bombles only
fit.bombus <- brm(bform, spec[spec$Genus == "Bombus",],
           cores=ncores,
           iter = 10^4,
           chains = 2,
           thin=1,
           inits=0,
           control = list(adapt_delta = 0.99))

write.ms.table(fit.bombus, "parasitism_bombus")

save(fit.bombus, spec,
     file="saved/parasiteBombusFitMod.Rdata")



## colonial only
fit.colonial <- brm(bform, spec[spec$Genus == "Bombus" |
                             spec$Genus == "Apis",],
           cores=ncores,
           iter = 10^4,
           chains = 2,
           thin=1,
           inits=0,
           control = list(adapt_delta = 0.99))

write.ms.table(fit.colonial, "parasitism_colonial")

save(fit.colonial, spec,
     file="saved/parasiteColonialFitMod.Rdata")
