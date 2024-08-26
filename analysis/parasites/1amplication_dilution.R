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

## site or lat as the geographic variable
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
vars_site <- c("Lat", "Area")
          
variables.to.log <- c("rare.degree", "Lat", "Net_BeeAbundance", "Area")

variables.to.log.1 <- c("Net_HBAbundance", "Net_BombusAbundance")

## loads specimen data
source("src/init.R")
## maybe remove SS? (only sampled one year). VC and UK were really
## grassy, odd meadows
spec.net <- filter(spec.net, Site != "VC" & Site != "UK")

spec.net$MeanITD[spec.net$Genus !=
                 "Bombus" & is.na(spec.net$Apidae)] <- NA
spec.net$rare.degree[!spec.net$Genus %in% c("Bombus", "Apis") &
                     is.na(spec.net$Apidae)] <- NA

## raw, non standardized data for plotting
spec.orig <- prepDataSEM(spec.net, variables.to.log, variables.to.log.1, 
                        standardize=FALSE)

## Make SEM weights and standardize data.
spec.net <- prepDataSEM(spec.net, variables.to.log, variables.to.log.1, 
                        vars_yearsr = vars_yearsr, vars_sp = vars_sp, 
                        vars_yearsrsp = vars_yearsrsp,
                        vars_site=vars_site)

##spec.net$Site <- factor(spec.net$Site, levels = c("JC", "SM", "SC", "MM", 
##                                                  "HM", "PL", "CH", "RP"))

spec.net$Site <- as.character(spec.net$Site)
## otherwise levels with no data are not properly dropped using subset
spec.net$Year <- as.character(spec.net$Year)
spec.net$GenusSpecies <- as.character(spec.net$GenusSpecies)

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

## define all the formulas for the different parts of the models
source("src/plant_poll_models.R")
## check ids
unique(spec.net$GenusSpecies[spec.net$Apidae == 1 &
                             is.na(spec.net$MeanITD)])
save(spec.net, spec.orig, file="saved/spec_weights.Rdata")

## Load phylogeny 
load("../../data/community_phylogeny.Rdata")
## Species that are not in the phylogeny are not used. brms is not
## allowing an incomplete phylogeny, to avoid the error we changed the
## species not present to one that is in the phylogeny.  We chose a
## species for which we did not do parasite screening and should not
## influence results.
not_in_phylo <- unique(spec.net$GenusSpecies[!spec.net$GenusSpecies
                                             %in%
                                             phylo$tip.label])
spec.bombus$GenusSpecies[spec.bombus$GenusSpecies %in%
                         not_in_phylo]<- "Agapostemon angelicus"

## **********************************************************
## Parasite models set up
## **********************************************************

## Multi species models
## social species abundances
xvars.multi.bombus.ss <-  c("Net_BeeDiversity",
                            "Net_BombusAbundance",
                            "Net_HBAbundance",
                            "rare.degree",
                            "MeanFloralDiversity",
                            "Lat",
                            "(1|Site)",
                            "(1|gr(GenusSpecies, cov = phylo_matrix))")

## bumble bee abundance only (vif sometimes indicate HB abundance and
## bombus abundance are colinear)
xvars.multi.bombus.ba <-  c("Net_BeeDiversity",
                            "Net_BombusAbundance",
                            "rare.degree",
                            "MeanFloralDiversity",
                            "Lat",
                            "(1|Site)",
                            "(1|gr(GenusSpecies, cov = phylo_matrix))")

## all bee abundance 
xvars.multi.bombus.all <-  c("Net_BeeDiversity",
                             "Net_BeeAbundance",
                             "rare.degree",
                             "MeanFloralDiversity",
                             "Lat",
                             "(1|Site)",
                             "(1|gr(GenusSpecies, cov = phylo_matrix))")

## single species models
## social species abundances
xvars.single.species.ss <-  c("Net_BeeDiversity",
                           "Net_BombusAbundance",
                           "Net_HBAbundance",
                           "rare.degree",
                           "MeanFloralDiversity",
                           "Lat",
                           "(1|Site)")

## all bee abundance 
xvars.single.species.all <-  c("Net_BeeDiversity",
                               "Net_BeeAbundance",
                               "rare.degree",
                               "MeanFloralDiversity",
                               "Lat",
                               "(1|Site)")

## **********************************************************
## community model, check assumptions first before adding parasites
## **********************************************************

run_plot_freq_model_diagnostics(remove_subset_formula(formula.flower.div),
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="students", site.lat=site.or.lat)

run_plot_freq_model_diagnostics(remove_subset_formula(formula.flower.abund),
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="students",site.lat=site.or.lat)

run_plot_freq_model_diagnostics(remove_subset_formula(formula.bee.abund),
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="students", site.lat=site.or.lat)

run_plot_freq_model_diagnostics(remove_subset_formula(formula.bombus.abund),
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="students", site.lat=site.or.lat)

run_plot_freq_model_diagnostics(remove_subset_formula(formula.HB.abund),
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="students", site.lat=site.or.lat)

run_plot_freq_model_diagnostics(remove_subset_formula(formula.bee.div),
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="gaussian",
                                site.lat=site.or.lat)

## **********************************************************
## Parasite presence
## **********************************************************
## bombus
table(spec.bombus$CrithidiaPresence[spec.bombus$WeightsPar == 1])
table(spec.bombus$ApicystisSpp[spec.bombus$WeightsPar ==1 ])

## social species abundance as x var
fit.bombus.ss <- runCombinedParasiteModels(spec.bombus,
                                        species.group="bombus",
                                        parasite =
                                            c("CrithidiaPresence",
                                              "ApicystisSpp"),
                                        xvars=xvars.multi.bombus.ss,
                                        ncores,
                                        iter = 10^4,
                                        chains = 1,
                                        thin=1,
                                        init=0, data2= list(phylo_matrix=phylo_matrix),
                                        SEM = TRUE, neg.binomial =  FALSE,
                                        site.lat=paste(site.or.lat,
                                                       "ss", sep="_"))

loo.crithidia.bombus.ss <- loo(fit.bombus.ss$fit, resp="CrithidiaPresence")
loo.apicystis.bombus.ss <- loo(fit.bombus.ss$fit, resp="ApicystisSpp")


## bombus abundance as x var
fit.bombus.ba <- runCombinedParasiteModels(spec.bombus,
                                        species.group="bombus",
                                        parasite =
                                            c("CrithidiaPresence",
                                              "ApicystisSpp"),
                                        xvars=xvars.multi.bombus.ba,
                                        ncores,
                                        iter = 10^4,
                                        chains = 1,
                                        thin=1,
                                        init=0, data2= list(phylo_matrix=phylo_matrix),
                                        SEM = TRUE, neg.binomial =  FALSE,
                                        site.lat=paste(site.or.lat,
                                                       "ba", sep="_"))

loo.crithidia.bombus.ba <- loo(fit.bombus.ba$fit, resp="CrithidiaPresence")
loo.apicystis.bombus.ba <- loo(fit.bombus.ba$fit, resp="ApicystisSpp")

## all bee abundance as xvar
fit.bombus.all <- runCombinedParasiteModels(spec.bombus,
                                        species.group="bombus",
                                        parasite =
                                            c("CrithidiaPresence",
                                              "ApicystisSpp"),
                                        xvars=xvars.multi.bombus.all,
                                        ncores,
                                        iter = 10^4,
                                        chains = 1,
                                        thin=1,
                                        init=0, data2= list(phylo_matrix=phylo_matrix),
                                        SEM = TRUE, neg.binomial =  FALSE,
                                        site.lat=paste(site.or.lat,
                                                       "all", sep="_"))

loo.crithidia.bombus.all <- loo(fit.bombus.all$fit, resp="CrithidiaPresence")
loo.apicystis.bombus.all <- loo(fit.bombus.all$fit,
                                resp="ApicystisSpp")

abundance.order <- c("social_species",
                     "bombus_abundance",
                     "all_bees")

loo.bombus.crithidia <- data.frame(Estimate=c(loo.crithidia.bombus.ss$estimates["looic", "Estimate"],
                                         loo.crithidia.bombus.ba$estimates["looic", "Estimate"],
                                         loo.crithidia.bombus.all$estimates["looic",
                                         "Estimate"]),
                                  SE=c(loo.crithidia.bombus.ss$estimates["looic", "SE"],
                                         loo.crithidia.bombus.ba$estimates["looic", "SE"],
                                         loo.crithidia.bombus.all$estimates["looic", "SE"]),
                                   abundVar= abundance.order)


## bombus and HB abundance has the best fit, but bombus and HB
## abundance are pretty colinear (VIF ~6). The next best fit is bombus
## abundance alone. 

loo.bombus.apicystis <- data.frame(Estimate=c(loo.apicystis.bombus.ss$estimates["looic", "Estimate"],
                                         loo.apicystis.bombus.ba$estimates["looic", "Estimate"],
                                         loo.apicystis.bombus.all$estimates["looic",
                                         "Estimate"]),
                                  SE=c(loo.apicystis.bombus.ss$estimates["looic", "SE"],
                                         loo.apicystis.bombus.ba$estimates["looic", "SE"],
                                         loo.apicystis.bombus.all$estimates["looic", "SE"]),
                                   abundVar= abundance.order)

## The best fit is the abundance of bombus and apis together, which in
## this model are not colinear. 

## **********************************************************
## Honey bees
table(spec.apis$CrithidiaPresence[spec.apis$WeightsPar == 1])
table(spec.apis$ApicystisSpp[spec.apis$WeightsPar ==1 ])

## social species abundance as x var
fit.apis.ss <- runCombinedParasiteModels(spec.apis, species.group="apis",
                                         parasite =
                                             c("CrithidiaPresence",
                                               "ApicystisSpp"),
                                         xvars=xvars.single.species.ss,
                                         ncores,
                                         iter = 10^4,
                                         chains = 1,
                                         thin=1,
                                         init=0,
                                         SEM = TRUE,
                                         neg.binomial = FALSE,
                                         site.lat=paste(site.or.lat,
                                                        "ss",
                                                        sep="_"))

loo.crithidia.apis.ss <- loo(fit.apis.ss$fit, resp="CrithidiaPresence")
loo.apicystis.apis.ss <- loo(fit.apis.ss$fit, resp="ApicystisSpp")

## all bee abundance as xvar
fit.apis.all <- runCombinedParasiteModels(spec.apis, species.group="apis",
                                          parasite =
                                              c("CrithidiaPresence",
                                                "ApicystisSpp"),
                                          xvars=xvars.single.species.all,
                                          ncores,
                                          iter = 10^4,
                                          chains = 1,
                                          thin=1,
                                          init=0,
                                          SEM = TRUE,
                                          neg.binomial = FALSE,
                                          site.lat=paste(site.or.lat,
                                                         "all",
                                                         sep="_"))


loo.crithidia.apis.all <- loo(fit.apis.all$fit, resp="CrithidiaPresence")
loo.apicystis.apis.all <- loo(fit.apis.all$fit, resp="ApicystisSpp")


abundance.order <- c("social_species",
                     "all_bees")

loo.apis.crithidia <- data.frame(Estimate=c(loo.crithidia.apis.ss$estimates["looic", "Estimate"],
                                         loo.crithidia.apis.all$estimates["looic",
                                         "Estimate"]),
                                  SE=c(loo.crithidia.apis.ss$estimates["looic", "SE"],
                                         loo.crithidia.apis.all$estimates["looic", "SE"]),
                                   abundVar= abundance.order)
loo.apis.crithidia

loo.apis.apicystis <- data.frame(Estimate=c(loo.apicystis.apis.ss$estimates["looic", "Estimate"],
                                         loo.apicystis.apis.all$estimates["looic",
                                         "Estimate"]),
                                  SE=c(loo.apicystis.apis.ss$estimates["looic", "SE"],
                                         loo.apicystis.apis.all$estimates["looic", "SE"]),
                                   abundVar= abundance.order)
loo.apis.apicystis

## models not distinguishable, equally good fit of bee abundance and
## HB + bombus abundance for both crithidia and apicystis 

## **********************************************************
## ## melissodes
## ## there are not enough positives to get this model to converge
## table(spec.melissodes$CrithidiaPresence[spec.melissodes$WeightsPar == 1])
## table(spec.melissodes$ApicystisSpp[spec.melissodes$WeightsPar ==1 ])

## fit.melissodes <- runCombinedParasiteModels(spec.melissodes,
##                                             species.group="melissodes",
##                                             parasite =
##                                                 c("CrithidiaPresence",
##                                                   "ApicystisSpp"),
##                                             xvars=xvars.single.species.all,
##                                             ncores,
##                                             iter = 10^4,
##                                             chains = 1,
##                                             thin=1,
##                                             init=0,
##                                             SEM = TRUE,
##                                             neg.binomial = FALSE,
##                                             site.lat=site.or.lat)

## ## anthophora
## ## there are not enough positives to get this model to converge
## table(spec.anthophora$CrithidiaPresence[spec.anthophora$WeightsPar == 1])
## table(spec.anthophora$ApicystisSpp[spec.anthophora$WeightsPar ==1 ])

## fit.anthophora <- runCombinedParasiteModels(spec.anthophora,
##                                             species.group="anthophora",
##                                             parasite =
##                                                 c("CrithidiaPresence",
##                                                   "ApicystisSpp"),
##                                             xvars=xvars.single.species.all,
##                                             ncores,
##                                             iter = 10^4,
##                                             chains = 1,
##                                             thin=1,
##                                             init=0,
##                                             SEM = TRUE,
##                                             neg.binomial = FALSE,
##                                             site.lat=site.or.lat)

## ## andrena
## ## there are not enough positives to get this model to converge
## table(spec.andrena$CrithidiaPresence[spec.andrena$WeightsPar == 1])
## table(spec.andrena$ApicystisSpp[spec.andrena$WeightsPar ==1 ])

## fit.andrena <- runCombinedParasiteModels(spec.andrena,
##                                             species.group="andrena",
##                                             parasite =
##                                                 c("CrithidiaPresence",
##                                                   "ApicystisSpp"),
##                                             xvars=xvars.single.species.all,
##                                             ncores,
##                                             iter = 10^4,
##                                             chains = 1,
##                                             thin=1,
##                                             init=0,
##                                             SEM = TRUE,
##                                             neg.binomial = FALSE,
##                                             site.lat=site.or.lat)
