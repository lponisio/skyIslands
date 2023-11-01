rm(list=ls())
## setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')
## setwd("C:/Users/na_ma/Dropbox (University of Oregon)/skyIslands")
setwd('~/Dropbox (University of Oregon)/skyislands')
ncores <- 1

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

## uses only net specimens, and drops syrphids
source("src/init.R")
## define all the formulas for the different parts of the models
source("src/plant_poll_models.R")
## check ids
unique(spec.net$GenusSpecies[spec.net$Apidae == 1 &
                             is.na(spec.net$MeanITD)])

## **********************************************************
## Parasite models set up
## **********************************************************
## Multi species models
xvars.multi.species <-  c("Net_NonBombusHBAbundance",
                          "Net_HBAbundance",
                          "Net_BombusAbundance",
                          "Net_BeeDiversity",
                          "rare.degree", "MeanITD",
                          "(1|Site)", "(1|GenusSpecies)")
## single species models
xvars.single.species <-  c("Net_NonBombusHBAbundance",
                           "Net_HBAbundance",
                           "Net_BombusAbundance",
                           "Net_BeeDiversity",
                           "rare.degree",
                           "(1|Site)")

## **********************************************************
## community model
## **********************************************************

bform.community <- bf.fabund + bf.fdiv +
  bf.babund 
  bf.bombusabund + bf.HBabund +
  bf.bdiv  +
  set_rescor(FALSE)

fit.community <- brm(bform.community, spec.net,
                     cores=ncores,
                     iter = 10^4,
                     chains = 1,
                     thin=1,
                     init=0,
                     control = list(adapt_delta = 0.99),
                     save_pars = save_pars(all = TRUE))
write.ms.table(fit.community,
               sprintf("parasitism_%s_%s",
                       species.group="all", parasite="none"))
r2 <- loo_R2(fit.community)
save(fit.community, spec.net, r2,
     file="saved/communityFit.Rdata")

## **********************************************************
## Parasite presence
## **********************************************************
## full model with all species, parasites
fit.all <- runCombinedParasiteModels(spec.all, species.group="all",
                                        parasites=c("CrithidiaPresence",
                                                    "ApicystisSpp"),
                                        xvars=xvars.multi.species,
                                        iter = 20^4,
                                        chains = 1,
                                        thin=1,
                                        init=0)

## bombus
fit.bombus <- runCombinedParasiteModels(spec.bombus, species.group="bombus",
                                        parasites=c("CrithidiaPresence",
                                                    "ApicystisSpp"),
                                        xvars=xvars.multi.species,
                                        iter = 10^4,
                                        chains = 1,
                                        thin=1,
                                        init=0)

## all apidae species (Bombus, Apis, Anthophora, Melissodes)
fit.apidae <- runCombinedParasiteModels(spec.apidae, species.group="apidae",
                                        parasites=c("CrithidiaPresence",
                                                    "ApicystisSpp"),
                                        xvars=xvars.multi.species,
                                        iter = 10^4,
                                        chains = 1,
                                        thin=1,
                                        init=0)

## Honey bees
fit.apis <- runCombinedParasiteModels(spec.apis, species.group="apis",
                                        parasites=c("CrithidiaPresence",
                                                    "ApicystisSpp"),
                                        xvars=xvars.single.species,
                                        iter = 10^4,
                                        chains = 1,
                                        thin=1,
                                        init=0)
