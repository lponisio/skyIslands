rm(list=ls())
## setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')
 setwd("C:/Users/na_ma/Dropbox (University of Oregon)/skyIslands")
## setwd('~/Dropbox (University of Oregon)/skyislands')
ncores <- 6

setwd("analysis/parasites")
source("src/misc.R")
source("src/writeResultsTable.R")
source("src/makeMultiLevelData.R")
source("src/runParasiteModels.R")
source("src/community_Model.R")
source("src/standardize_weights.R")

## all of the variables that are explanatory variables and thus need
## to be centered
vars_yearsr <- c("MeanFloralAbundance",
          "MeanFloralDiversity",
          "Net_BeeDiversity",
          "Lat", "SRDoy"  
          )
vars_yearsrsp <- "rare.degree"
vars_sp <- "MeanITD"
          

variables.to.log <- c("rare.degree", "MeanITD")

variables.to.log.1<- c("Net_HBAbundance", "Net_BombusAbundance", 
                       "Net_NonBombusHBAbundance")
                      

## uses only net specimens, and drops syrphids
source("src/init.R")
## Make SEM weights and standardize data.
spec.net <- prepDataSEM(spec.net, variables.to.log, variables.to.log.1, 
                        vars_yearsr = vars_yearsr, vars_sp = vars_sp, 
                        vars_yearsrsp = vars_yearsrsp)

## bombus only data
spec.bombus <- spec.net
spec.bombus$WeightsPar[spec.bombus$Genus != "Bombus"] <- 0

## apis only data
spec.apis <- spec.net
spec.apis$WeightsPar[spec.apis$Genus != "Apis"] <- 0

## melissodes only data
spec.melissodes <- spec.net
spec.melissodes$WeightsPar[spec.melissodes$Genus != "Melissodes"] <- 0

spec.apidae <- spec.net
spec.apidae$WeightsPar[spec.apidae$Family != "Apidae"] <- 0
## define all the formulas for the different parts of the models
source("src/plant_poll_models.R")
## check ids
unique(spec.net$GenusSpecies[spec.net$Apidae == 1 &
                             is.na(spec.net$MeanITD)])

## Load phylogeny 
load("../../data/community_phylogeny.Rdata")
## Species that are not in the phylogeny are not used. brms is not allowing an incomplete
## phylogeny, to avoid the error we changed the species not present to one that is in the phylogeny. 
## We chose a species for which we did not do parasite screening and should not influence results.
not_in_phylo <- unique(spec.net$GenusSpecies[!spec.net$GenusSpecies %in% phylo$tip.label])
spec.bombus$GenusSpecies[spec.bombus$GenusSpecies %in% not_in_phylo]<- "Agapostemon angelicus"

## **********************************************************
## Parasite models set up
## **********************************************************
## Multi species models
xvars.multi.bombus <-  c("Net_BombusAbundance",
                          "Net_BeeDiversity",
                          "rare.degree", "MeanITD", "(1|Site)",
                           "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.multi.species <-  c("Net_NonBombusHBAbundance",
                          "Net_BeeDiversity",
                          "rare.degree", "MeanITD",
                          "(1|Site)", 
                           "(1|GenusSpecies)")


## single species models
xvars.single.species <-  c("Net_HBAbundance",
                           "Net_BeeDiversity",
                           "rare.degree",
                           "(1|Site)")

## **********************************************************
## community model
## **********************************************************

bform.community <- bf.fabund + bf.fdiv +
  bf.babund +
  bf.bombusabund + bf.HBabund +
  bf.bdiv  +
  set_rescor(FALSE)

fit.community <- runCommunityModels(bform.community, spec.net, 
                   ncores, 
                   iter = 10^4,
                   chains = 1,
                   thin=1,
                   init=0)



## **********************************************************
## Parasite presence
## **********************************************************
## full model with Melissodes, parasites


fit.parasites <- runCombinedParasiteModels(spec.net, species.group="melissodes",
                             parasite= c("CrithidiaPresence",
                                         "ApicystisSpp"),
                             xvars=xvars.multi.species,
                             iter = 10^4,
                             chains = 1,
                             thin=1,
                             init=0)

## bombus

fit.bombus <- runCombinedParasiteModels(
  spec.bombus, species.group="bombus",
                                        parasite = c("CrithidiaPresence", "ApicystisSpp"),
                                        xvars=xvars.multi.bombus,
                                        ncores,
                                        iter = 10^4,
                                        chains = 1,
                                        thin=1,
                                        init=0, data2= list(phylo_matrix=phylo_matrix))



## all apidae species (Bombus, Apis, Anthophora, Melissodes)
fit.apidae <- runCombinedParasiteModels(spec.apidae, species.group="apidae",
                                        parasites=c("CrithidiaPresence",
                                                    "ApicystisSpp"),
                                        xvars=xvars.multi.species,
                                        iter = 2*(10^4),
                                        chains = 1,
                                        thin=1,
                                        init=0)



## Honey bees
fit.apis <- runCombinedParasiteModels(spec.apis, species.group="apis",
                                        parasites=c("CrithidiaPresence",
                                                    "ApicystisSpp"),
                                        xvars=xvars.single.species,
                                        iter = 2*(10^4),
                                        chains = 1,
                                        thin=1,
                                        init=0)


