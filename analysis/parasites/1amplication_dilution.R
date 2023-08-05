rm(list=ls())
## setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')
 setwd("C:/Users/na_ma/Dropbox (University of Oregon)/skyIslands")
##setwd('~/Dropbox (University of Oregon)/skyislands')
ncores <- 2

setwd("analysis/parasites")
source("src/misc.R")
source("src/writeResultsTable.R")
source("src/makeMultiLevelData.R")
source("src/runParasiteModels.R")
## all of the variables that are explanatory variables and thus need
## to be centered
vars <- c("MeanFloralAbundance",
          "MeanFloralDiversity",
          "Net_BeeDiversity",
          "Lat", "SRDoy",
           "MeanITD",
          ## "r.degree", ## across all networks
          "rare.degree"  ## site-year level
          )
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
                          ## "Net_HBAbundance",
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
## Parasite presence
## **********************************************************
## full model with all species, parasites

 fit <- runParasiteModels(spec.data= spec.all,
                         species.group="all",
                         parasite="ParasitePresence",
                          xvars= xvars.multi.species)

## **********************************************************
## Crithidia models
## **********************************************************

## ## bombus only models
## fit.bombus.cbombi <- runParasiteModels(spec.bombus, "bombus",
##                                        "CrithidiaBombi",
##                                        xvars.multi.species)
## fit.bombus.cexpoeki <- runParasiteModels(spec.bombus, "bombus",
##                                          "CrithidiaExpoeki",
##                                          xvars.multi.species)
 fit.bombus.cspp <- runParasiteModels(spec.bombus, "bombus", "CrithidiaPresence",
                                      xvars.multi.species)

## ## all apidae species (Bombus, Apis, Anthophora, Melissodes)
## fit.apidae.cbombi <- runParasiteModels(spec.apidae, "apidae",
##                                        "CrithidiaBombi",
##                                        xvars.multi.species)
## fit.apidae.cexpoeki <- runParasiteModels(spec.apidae, "apidae",
##                                          "CrithidiaExpoeki",
##                                          xvars.multi.species)
## fit.apidae.cspp <- runParasiteModels(spec.apidae, "apidae",
##                                      "CrithidiaPresence",
##                                      xvars.multi.species)

## ## Apis only
## fit.apis.cbombi <- runParasiteModels(spec.apis, "apis",
##                                      "CrithidiaBombi",
##                                      xvars.single.species)
## fit.apis.cexpoeki <- runParasiteModels(spec.apis, "apis",
##                                        "CrithidiaExpoeki",
##                                        xvars.single.species)
## fit.apis.cspp <- runParasiteModels(spec.apis, "apis",
##                                    "CrithidiaPresence",
##                                    xvars.single.species)

## **********************************************************
## ApicystisSpp models
## **********************************************************
## bombus only models
fit.bombus.apicystis <- runParasiteModels(spec.bombus, "bombus",
                                          "ApicystisSpp",
                                          xvars.multi.species)
## all apidae species (Bombus, Apis, Anthophora, Melissodes)
fit.apidae.apicystis <- runParasiteModels(spec.apidae, "apidae",
                                          "ApicystisSpp",
                                          xvars.multi.species)
## Apis only
fit.apis.apicystis <- runParasiteModels(spec.apis, "apis",
                                        "ApicystisSpp",
                                        xvars.single.species)
