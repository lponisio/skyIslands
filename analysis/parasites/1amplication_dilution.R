rm(list=ls())
## setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')
 setwd("C:/Users/na_ma/Dropbox (University of Oregon)/skyIslands")
## setwd('~/Dropbox (University of Oregon)/skyislands')
ncores <- 5

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
          "Year", "SRDoy"  
          )
vars_sp <- c("MeanITD",
          "rare.degree")

variables.to.log <- c("rare.degree", "MeanITD")

variables.to.log.1<- c("Net_HBAbundance", "Net_BombusAbundance", 
                       "Net_NonBombusHBAbundance")
                      

## uses only net specimens, and drops syrphids
source("src/init.R")
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
                          "rare.degree", "MeanITD",
                          "(1|Site)", "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.multi.species <-  c("Net_NonBombusHBAbundance",
                          "Net_BeeDiversity",
                          "rare.degree", "MeanITD",
                          "(1|Site)", "(1|GenusSpecies)")


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

fit.community <- brm(bform.community, spec.net,
                     cores=ncores,
                     iter = 10^4,
                     chains = 1,
                     thin=1,
                     control = list(adapt_delta = 0.99),
                     save_pars = save_pars(all = TRUE))
write.ms.table(fit.community,
               sprintf("parasitism_%s_%s",
                       species.group="all", parasite="none"))
r2loo <- loo_R2(fit.community)
r2 <- rstantools::bayes_R2(fit.community)
save(fit.community, spec.net, r2, spec.orig,
     file="saved/communityFit.Rdata")

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

fit.bombus <- runCombinedParasiteModels(spec.bombus, species.group="bombus",
                                        parasite = c("CrithidiaPresence", "ApicystisSpp"),
                                        xvars=xvars.multi.bombus,
                                        ncores,
                                        iter = 2*(10^4),
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


