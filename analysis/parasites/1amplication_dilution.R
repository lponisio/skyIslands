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
source("src/community_Model.R")
source("src/standardize_weights.R")
source("src/runPlotFreqModelDiagnostics.R")

site.or.lat <- "lat"

## all of the variables that are explanatory variables and thus need
## to be centered
vars_yearsr <- c("MeanFloralAbundance",
                 "MeanFloralDiversity",
                 "Net_BeeDiversity",
                 "Net_BeeAbundance",
                 "Net_BombusAbundance",
                 "Net_HBAbundance",
                 "Net_NonBombusHBAbundance",
                 "SRDoyPoly1",
                 "SRDoyPoly2"
                 )
vars_yearsrsp <- "rare.degree"
vars_sp <- "MeanITD"
vars_site <- "Lat"
          
variables.to.log <- c("rare.degree", "MeanITD", "Lat",
                      "Net_NonBombusHBAbundance", "Net_BeeAbundance")

variables.to.log.1 <- c("Net_HBAbundance", "Net_BombusAbundance")

## loads specimen data
source("src/init.R")
spec.net <- filter(spec.net, Site != "VC" & Site != "UK" & Site != "SS")


## Make SEM weights and standardize data.
spec.net <- prepDataSEM(spec.net, variables.to.log, variables.to.log.1, 
                        vars_yearsr = vars_yearsr, vars_sp = vars_sp, 
                        vars_yearsrsp = vars_yearsrsp,
                        vars_site=vars_site)
spec.net$Site <- factor(spec.net$Site, levels = c("JC", "SM", "SC", "MM", 
                                                  "HM", "PL", "CH", "RP"))
spec.net$Year <- as.factor(spec.net$Year)  
## bombus only data
spec.bombus <- spec.net
spec.bombus$WeightsPar[spec.bombus$Genus != "Bombus"] <- 0
spec.bombus$WeightsSp[spec.bombus$Genus != "Bombus"] <- 0

## apis only data
spec.apis <- spec.net
spec.apis$WeightsPar[spec.apis$Genus != "Apis"] <- 0
spec.apis$WeightsSp[spec.apis$Genus != "Apis"] <- 0
## melissodes only data
spec.melissodes <- spec.net
spec.melissodes$WeightsPar[spec.melissodes$Genus != "Melissodes"] <- 0
spec.melissodes$WeightsSp[spec.melissodes$Genus != "Melissodes"] <- 0

spec.apidae <- spec.net
spec.apidae$WeightsPar[spec.apidae$Family != "Apidae"] <- 0
spec.apidae$WeightsSp[spec.apidae$Family != "Apidae"] <- 0
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
xvars.multi.bombus <-  c("Net_BeeDiversity", "Net_BombusAbundance",
                          "rare.degree",  "(1|Site)",
                           "(1|gr(GenusSpecies, cov = phylo_matrix))")


## For apicystis MeanITD makes sense but for crithidia it decreases the fit
xvars.multi.bombus2 <-  c("Net_BeeDiversity", "Year",
                         "rare.degree", "MeanITD", "(1|Site)",
                         "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.multi.species <-  c("Year",
                          "Net_BeeDiversity",
                          "Net_BeeAbundance",
                          "rare.degree",
                          "(1|Site)", 
                           "(1|GenusSpecies)")

## single species models
xvars.single.species <-  c("Year",
                           "Net_BeeDiversity",
                           "Net_BeeAbundance",
                           "rare.degree",
                           "(1|Site)")

## Apis
xvars.single.species <-  c("Year",
                           "Net_BeeDiversity",
                           "Net_HBAbundance",
                           "rare.degree",
                           "(1|Site)")


## **********************************************************
## community model, check assumptions first before adding parasites
## **********************************************************

run_plot_freq_model_diagnostics(formula.flower.div.nosub,
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="students")

run_plot_freq_model_diagnostics(formula.flower.abund.nosub,
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="students")

run_plot_freq_model_diagnostics(formula.bee.abund.nosub,
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="students")

run_plot_freq_model_diagnostics(formula.bombus.abund.nosub,
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="students")

run_plot_freq_model_diagnostics(formula.HB.abund.nosub,
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="students")

run_plot_freq_model_diagnostics(formula.bee.div.nosub,
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="gaussian")

## **********************************************************
## Parasite presence
## **********************************************************
## full model with Melissodes 
fit.melissodes <- runCombinedParasiteModels(spec.melissodes, species.group="melissodes",
                                            parasite = c("CrithidiaPresence", "ApicystisSpp"),
                                            xvars=xvars.single.species,
                                            ncores,
                                            iter = 10^4,
                                            chains = 1,
                                            thin=1,
                                            init=0,
                                            SEM = TRUE, neg.binomial =
                                                            TRUE,
                                            site.lat=site.or.lat)


## bombus
fit.bombus <- runCombinedParasiteModels(spec.bombus, species.group="bombus",
                                        parasite = c("CrithidiaPresence", "ApicystisSpp"),
                                        xvars=xvars.multi.bombus,
                                        ncores,
                                        iter = 10^4,
                                        chains = 1,
                                        thin=1,
                                        init=0, data2= list(phylo_matrix=phylo_matrix),
                                        SEM = TRUE, neg.binomial =
                                                        TRUE,
                                        site.lat=site.or.lat)


## Honey bees
fit.apis <- runCombinedParasiteModels(spec.melissdoes, species.group="apis",
                                            parasite = c("CrithidiaPresence", "ApicystisSpp"),
                                            xvars=xvars.single.species,
                                            ncores,
                                            iter = 10^4,
                                            chains = 1,
                                            thin=1,
                                            init=0,
                                            SEM = TRUE, neg.binomial =
                                                            TRUE,
                                            site.lat=site.or.lat)


