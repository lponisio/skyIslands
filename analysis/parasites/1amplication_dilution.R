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
## maybe remove SS? (only sampled one year). VC and UK were really
## grassy, odd meadows

## spec.net <- filter(spec.net, Site != "VC" & Site != "UK" & Site != "SS")
spec.net <- filter(spec.net, Site != "VC")


## because only Bombus and apis models converge, setted the rest of
## the trait data to NA so that the variables scale properly
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

## can repeat for other genera but have deleted since the models do
## not converge

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
xvars.ba <-  c("Net_BeeDiversity",
                            "Net_BombusAbundance",
                            "rare.degree",
                            "MeanFloralDiversity",
                            "Lat",
                            "(1|Site)",
                            "(1|gr(GenusSpecies, cov = phylo_matrix))")


## honey bee abundance only (vif sometimes indicates floral diversity and
## bombus abundance are colinear)
xvars.ha <-  c("Net_BeeDiversity",
                            "Net_HBAbundance",
                            "rare.degree",
                            "MeanFloralDiversity",
                            "Lat",
                            "(1|Site)",
                            "(1|gr(GenusSpecies, cov = phylo_matrix))")



## all bee abundance 
xvars.all <-  c("Net_BeeDiversity",
                             "Net_BeeAbundance",
                             "rare.degree",
                             "MeanFloralDiversity",
                             "Lat",
                             "(1|Site)",
                             "(1|gr(GenusSpecies, cov = phylo_matrix))")


## diversity only 
xvars.div <-  c("Net_BeeDiversity",
                "rare.degree",
                "MeanFloralDiversity",
                "Lat",
                "(1|Site)",
                "(1|gr(GenusSpecies, cov = phylo_matrix))")



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



runCombinedParasiteModesWrapper <- function(data,
                                            nvar.name,
                                            xvars,
                                            species.group=species.group,
                                            iste.lat=site.lat,
                                            ncores=ncores,
                                            iter=10^4,
                                            chains=1,
                                            SEM = TRUE,
                                            neg.binomial =  FALSE,
                                            data2= list(phylo_matrix=phylo_matrix),
                                            ){
    fit <- runCombinedParasiteModels(spec.bombus,
                                     species.group="bombus",
                                     parasite =
                                         c("CrithidiaPresence",
                                           "ApicystisSpp"),
                                     xvars=xvars,
                                     ncores,
                                     iter = iter,
                                     chains = chains,
                                     thin=1,
                                     init=0,
                                     data2= data2,
                                     SEM = SEM,
                                     neg.binomial =  neg.binomial,
                                     site.lat=paste(site.or.lat,
                                                    xvar.name, sep="_"))
    
    loo.crithidia <- loo(fit$fit, resp="CrithidiaPresence")
    loo.apicystis <- loo(fit$fit, resp="ApicystisSpp")

    return(list(fit=fit, loo.crithidia= loo.crithidia,
                loo.apicystis=loo.apicystis))

}

## moel running function also runs frequentist models and
## diagnostics. Bayesian models are beta binomial, frequentist are
## binomial because there are no packages that use beta binomial

## social species abundance as x var
fit.bombus.ss <- runCombinedParasiteModels(spec.bombus,
                                        species.group="bombus",
                                        parasite =
                                            c("CrithidiaPresence",
                                              "ApicystisSpp"),
                                        xvars=xvars.ss,
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
                                        xvars=xvars.ba,
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
                                        xvars=xvars.all,
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


loo.bombus.crithidia <- makeLooTable(parasite="CrithidiaSpp",
                                     genus="Bombus",
                                     abundance.order=abundance.order,
                                     list(loo.crithidia.bombus.ss,
                                          loo.crithidia.bombus.ba,
                                          loo.crithidia.bombus.all)
                                     )

## not including ss because HB abundance and bombus abundance are colinear
loo_compare(loo.crithidia.bombus.ba, loo.crithidia.bombus.all)

## bombus and HB abundance has the best fit, but bombus and HB
## abundance are pretty colinear (VIF ~6). The next best fit is bombus
## abundance alone. 

loo.bombus.apicystis <- makeLooTable(parasite="ApicystisSpp",
                                     genus="Bombus",
                                     abundance.order=abundance.order,
                                     list(loo.apicystis.bombus.ss,
                                          loo.apicystis.bombus.ba,
                                          loo.apicystis.bombus.all)
                                     )

loo_compare(loo.apicystis.bombus.ss, loo.apicystis.bombus.ba,
            loo.apicystis.bombus.all)

## The best fit is the abundance of bombus and apis together, which in
## this model are not colinear. 

## **********************************************************
## Honey bees
## **********************************************************
table(spec.apis$CrithidiaPresence[spec.apis$WeightsPar == 1])
table(spec.apis$ApicystisSpp[spec.apis$WeightsPar ==1 ])

## social species abundance as x var
fit.apis.ss <- runCombinedParasiteModels(spec.apis, species.group="apis",
                                         parasite =
                                             c("CrithidiaPresence",
                                               "ApicystisSpp"),
                                         xvars=xvars.ss[-length(xvars.ss)],
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

## HB  abundance as x var
fit.apis.ha <- runCombinedParasiteModels(spec.apis, species.group="apis",
                                         parasite =
                                             c("CrithidiaPresence",
                                               "ApicystisSpp"),
                                         xvars=xvars.ha[-length(xvars.ha)],
                                         ncores,
                                         iter = 10^4,
                                         chains = 1,
                                         thin=1,
                                         init=0,
                                         SEM = TRUE,
                                         neg.binomial = FALSE,
                                         site.lat=paste(site.or.lat,
                                                        "ha",
                                                        sep="_"))

loo.crithidia.apis.ha <- loo(fit.apis.ha$fit, resp="CrithidiaPresence")
loo.apicystis.apis.ha <- loo(fit.apis.ha$fit, resp="ApicystisSpp")

## all bee abundance as xvar
fit.apis.all <- runCombinedParasiteModels(spec.apis, species.group="apis",
                                          parasite =
                                              c("CrithidiaPresence",
                                                "ApicystisSpp"),
                                          xvars=xvars.all[-length(xvars.all)],
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



## no bee abundance only bee diversity  as xvar
fit.apis.div <- runCombinedParasiteModels(spec.apis, species.group="apis",
                                          parasite =
                                              c("CrithidiaPresence",
                                                "ApicystisSpp"),
                                          xvars=xvars.div[-length(xvars.div)],
                                          ncores,
                                          iter = 10^4,
                                          chains = 1,
                                          thin=1,
                                          init=0,
                                          SEM = TRUE,
                                          neg.binomial = FALSE,
                                          site.lat=paste(site.or.lat,
                                                         "div",
                                                         sep="_"))


loo.crithidia.apis.div <- loo(fit.apis.div$fit, resp="CrithidiaPresence")
loo.apicystis.apis.div <- loo(fit.apis.div$fit, resp="ApicystisSpp")


## loo combined
abundance.order <- c("social_species",
                     "HB_abundance",
                     "all_bees",
                     "diversity_noabund")

loo.apis.crithidia <- makeLooTable(parasite="CrithidiaSpp",
                                     genus="Apis",
                                     abundance.order=abundance.order,
                                     list(loo.crithidia.apis.ss,
                                          loo.crithidia.apis.ha,
                                          loo.crithidia.apis.all,
                                          loo.crithidia.apis.div)
                                     )

loo_compare(loo.crithidia.apis.ss, loo.crithidia.apis.ha,
            loo.crithidia.apis.all, loo.crithidia.apis.div)

abundance.order <- c("social_species",
                     "all_bees")
loo.apis.apicystis <- makeLooTable(parasite="ApicystisSpp",
                                     genus="Apis",
                                     abundance.order=abundance.order,
                                     list(loo.apicystis.apis.ss,
                                          loo.apicystis.apis.all)
                                     )

loo_compare(loo.apicystis.apis.ss, loo.apicystis.apis.all)

## models not distinguishable, equally good fit of bee abundance and
## HB + bombus abundance for both crithidia and apicystis 

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
##                                             xvars=xvars.ss[-length(xvars.ss)],
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
##                                             xvars=xvars.ss[-length(xvars.ss)],
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
##                                             xvars=xvars.ss[-length(xvars.ss)],
##                                             ncores,
##                                             iter = 10^4,
##                                             chains = 1,
##                                             thin=1,
##                                             init=0,
##                                             SEM = TRUE,
##                                             neg.binomial = FALSE,
##                                             site.lat=site.or.lat)
