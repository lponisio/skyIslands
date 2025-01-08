rm(list=ls())
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
                 "SRDoyPoly1",
                 "SRDoyPoly2"
                 )
vars_yearsrsp <- "rare.degree"
vars_sp <- "MeanITD"
vars_site <- c("Lat", "Area")

variables.to.log <- c("rare.degree", "Lat", "Net_BeeAbundance", "Area")

## some zeros in data
variables.to.log.1 <- c("Net_HBAbundance", "Net_BombusAbundance")

## loads specimen data
source("src/init.R")

## drop VC (Valles caldera) because it was more of a grassland, we only surveyed it one year)
print("Before dropping VC")
dim(spec.net)
spec.net <- filter(spec.net, Site != "VC")
print("After dropping VC")
dim(spec.net)
## because only Bombus and apis models converge, setted the rest of
## the trait data to NA so that the variables scale properly
screened.bombus <- unique(spec.net$GenusSpecies[spec.net$Apidae == 1 &
                                                spec.net$Genus == "Bombus"])
screened.bombus <- screened.bombus[!is.na(screened.bombus)]

spec.net$MeanITD[!spec.net$GenusSpecies %in% screened.bombus] <- NA

spec.net$rare.degree[!spec.net$GenusSpecies %in%
                     c("Apis mellifera", screened.bombus)] <- NA

dim(spec.net)
## raw, non standardized data for plotting
spec.orig <- prepDataSEM(spec.net, variables.to.log, variables.to.log.1,
                         standardize=FALSE)

## Make SEM weights and standardize data.
spec.net <- prepDataSEM(spec.net, variables.to.log, variables.to.log.1,
                        vars_yearsr = vars_yearsr, vars_sp = vars_sp,
                        vars_yearsrsp = vars_yearsrsp,
                        vars_site=vars_site)

spec.net$Site <- as.character(spec.net$Site)
## otherwise levels with no data are not properly dropped using subset
spec.net$Year <- as.character(spec.net$Year)
spec.net$GenusSpecies <- as.character(spec.net$GenusSpecies)

## bombus only data
spec.bombus <- makeGenusSubset(spec.net, "Bombus")
## apis only data
spec.apis <- makeGenusSubset(spec.net, "Apis")
## can repeat for other genera but have deleted since the models do
## not converge

## define all the formulas for the different parts of the models
source("src/plant_poll_models.R")

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

## only bombus model is muti species given the others do not converge
spec.bombus$GenusSpecies[spec.bombus$GenusSpecies %in%
                         not_in_phylo]<- "Agapostemon angelicus"

## **********************************************************
## Parasite models set up
## **********************************************************
## phylogeny must be last in all xvar sets
## social species abundances
xvars.ss <-  c("Net_BeeDiversity",
               "Net_BombusAbundance",
               "Net_HBAbundance",
               "rare.degree",
               "MeanFloralDiversity",
               "Lat",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

## bumble bee abundance only (vif sometimes indicates HB abundance and
## bombus abundance are colinear)
xvars.ba <- xvars.ss[xvars.ss != "Net_HBAbundance"]

## honey bee abundance only (vif sometimes indicates floral diversity and
## bombus abundance are colinear)
xvars.ha <- xvars.ss[xvars.ss != "Net_BombusAbundance"]

## all bee abundance
xvars.all <- c("Net_BeeAbundance",
               xvars.ss[!xvars.ss %in%
                        c("Net_HBAbundance", "Net_BombusAbundance")])

## diversity only
xvars.div <- xvars.all[xvars.all != "Net_BeeAbundance"]

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
## Bombus
## **********************************************************
table(spec.bombus$CrithidiaPresence[spec.bombus$WeightsPar == 1])
table(spec.bombus$ApicystisSpp[spec.bombus$WeightsPar ==1 ])

## model running function also runs frequentist models and
## diagnostics. Bayesian models are beta binomial, frequentist are
## binomial because there are no packages that include beta binomial

xvar.order <- c("social_species",
                "bombus_abundance",
                "all_bees")

## social species abundance as x var
bombus.ss <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.ss,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       site.lat=site.or.lat,
                                       xvar.name=xvar.order[1])
## HB and bombus abundance are colinear in crithidia model, so not a valid model

## bombus abundance only
bombus.ba <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.ba,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       site.lat=site.or.lat,
                                       xvar.name=xvar.order[2])

## all bee abundance
bombus.all <- runCombinedParasiteModels(spec.data= spec.bombus,
                                        species.group="bombus",
                                        xvars=xvars.all,
                                        ncores=ncores,
                                        data2= list(phylo_matrix=phylo_matrix),
                                        site.lat=site.or.lat,
                                        xvar.name=xvar.order[3])

## **********************************************************
## Bombus loo summaries
## **********************************************************
## crithidia
## **********************************************************
## not including ss because HB abundance and bombus abundance are colinear
load("saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat_all_bees.Rdata")
bombus.all <- fit.parasite
load("saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat_bombus_abundance.Rdata")
bombus.ba <- fit.parasite

loo.crithidia.all <- loo(bombus.all, resp="CrithidiaPresence")
loo.crithidia.ba <- loo(bombus.ba, resp="CrithidiaPresence")

# bombus.loo.crithidia <- list(ba=bombus.ba$loo.crithidia,
#                              all=bombus.all$loo.crithidia)
# 
# sum.loo.bombus.crithidia <- makeLooTable(parasite="CrithidiaSpp",
#                                          genus="Bombus",
#                                          abundance.order=xvar.order[-1],
#                                          bombus.loo.crithidia
#                                          )
loo_compare(loo.crithidia.all, loo.crithidia.ba)
## bombus abundance an all abundance model fits are not distinguishable

## **********************************************************
## apicystis
## **********************************************************
## not including ba (bombus abundance) since HB and abundance
## abundance are not colinear in these models
load("saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat_social_species.Rdata")
bombus.ss <- fit.parasite

loo.apicystis.all <- loo(bombus.all, resp="ApicystisSpp")
loo.apicystis.ss <- loo(bombus.ss, resp="ApicystisSpp")

# bombus.loo.apicystis <- list(ss=bombus.ss$loo.apicystis,
#                              all=bombus.all$loo.apicystis)
# 
# sum.loo.bombus.apicystis <- makeLooTable(parasite="ApicystisSpp",
#                                          genus="Bombus",
#                                          abundance.order=xvar.order[-2],
#                                          bombus.loo.apicystis
#                                          )
loo_compare(loo.apicystis.all, loo.apicystis.ss)
## The best fit is the abundance of bombus and apis together, which in
## this model are not colinear.

## **********************************************************
## Honey bees
## **********************************************************
table(spec.apis$CrithidiaPresence[spec.apis$WeightsPar == 1])
table(spec.apis$ApicystisSpp[spec.apis$WeightsPar ==1 ])


xvar.order <- c("social_species",
                "apis_abundance",
                "all_bees",
                "diversity")

## social species abundance as x var
apis.ss <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.ss[-length(xvars.ss)],
                                     ncores=ncores,
                                     site.lat=site.or.lat,
                                     xvar.name=xvar.order[1])
## HB and apis abundance are colinear in crithidia model, so not a valid model

## apis abundance only
apis.ha <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.ha[-length(xvars.ha)],
                                     ncores=ncores,
                                     site.lat=site.or.lat,
                                     xvar.name=xvar.order[2])

## all bee abundance
apis.all <- runCombinedParasiteModels(spec.data= spec.apis,
                                      species.group="apis",
                                      xvars=xvars.all[-length(xvars.all)],
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      site.lat=site.or.lat,
                                      xvar.name=xvar.order[3])

## no abundances, just bee diversity (abundances are colinear in some models)
apis.div <- runCombinedParasiteModels(spec.data= spec.apis,
                                      species.group="apis",
                                      xvars=xvars.div[-length(xvars.div)],
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      site.lat=site.or.lat,
                                      xvar.name=xvar.order[4])

## **********************************************************
## Apis loo summaries
## **********************************************************
## crithidia
## **********************************************************
## not including ss because of colinearity
load("saved/parasiteFit_apis_CrithidiaPresenceApicystisSpp_lat_all_bees.Rdata")
apis.all <- fit.parasite
load("saved/parasiteFit_apis_CrithidiaPresenceApicystisSpp_lat_apis_abundance.Rdata")
apis.ha <- fit.parasite
load("saved/parasiteFit_apis_CrithidiaPresenceApicystisSpp_lat_diversity.Rdata")
apis.div <- fit.parasite

# apis.loo.crithidia <- list(ha=apis.ha$loo.crithidia,
#                            all=apis.all$loo.crithidia,
#                            div=apis.div$loo.crithidia)
# 
# sum.loo.apis.crithidia <- makeLooTable(parasite="CrithidiaSpp",
#                                        genus="Apis",
#                                        abundance.order=xvar.order[-1],
#                                        apis.loo.crithidia
#                                        )
loo.crithidia.all <- loo(apis.all, resp="CrithidiaPresence")
loo.crithidia.ha <- loo(apis.ha, resp="CrithidiaPresence")
loo.crithidia.div <- loo(apis.div, resp="CrithidiaPresence")

loo_compare(loo.crithidia.all, loo.crithidia.ha, loo.crithidia.div)

## **********************************************************
## apicystis
## **********************************************************
# apis.loo.apicystis <- list(ha=apis.ha$loo.apicystis,
#                            all=apis.all$loo.apicystis,
#                            div=apis.div$loo.apicystis)
# 
# sum.loo.apis.apicystis <- makeLooTable(parasite="ApicystisSpp",
#                                        genus="Apis",
#                                        abundance.order=xvar.order[-1],
#                                        apis.loo.apicystis
#                                        )

loo.apicystis.all <- loo(apis.all, resp="ApicystisSpp")
loo.apicystis.ha <- loo(apis.ha, resp="ApicystisSpp")
loo.apicystis.div <- loo(apis.div, resp="ApicystisSpp")

loo_compare(loo.apicystis.all, loo.apicystis.ha, loo.apicystis.div)

## **********************************************************
## save models and loo results
## **********************************************************

write.csv(rbind(loo.bombus.crithidia,
                loo.apis.crithidia,
                loo.bombus.apicystis,
                loo.apis.apicystis),
          row.names=FALSE,
          file="saved/loo.csv")

write.table(rbind(loo.bombus.crithidia,
                  loo.apis.crithidia,
                  loo.bombus.apicystis,
                  loo.apis.apicystis),
            row.names=FALSE,
            file="saved/loo.txt", sep= " & ")

save(loo.apis.apicystis,
     loo.apis.crithidia,
     loo.bombus.apicystis,
     loo.bombus.crithidia,
     fit.apis.all,
     fit.bombus.all,
     fit.apis.ss,
     fit.bombus.ss,
     fit.apis.ha,
     fit.bombus.ba,
     file="saved/all_models_loo.R")
