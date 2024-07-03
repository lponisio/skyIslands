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
# spec.net$Site <- factor(spec.net$Site, levels = c("JC", "SM", "SC", "MM", 
                                                  "HM", "PL", "CH", "RP"))

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
xvars.multi.bombus <-  c("Net_BeeDiversity*Year", "Net_BombusAbundance",
                         "rare.degree",
                         "MeanFloralAbundance",
                         "MeanFloralDiversity",
                         "MeanITD",
                         "(1|Site)",
                         "(1|gr(GenusSpecies, cov = phylo_matrix))")

## mean ITD causing problems

## single species models
xvars.single.species <-  c("Net_BeeDiversity",
                           "Net_BeeAbundance",
                           "MeanFloralAbundance",
                           "MeanFloralDiversity",
                           "rare.degree",
                           "(1|Site)")

## Apis
xvars.apis <-  c("Net_BeeDiversity",
                 "Net_HBAbundance",
                 "MeanFloralAbundance",
                 "MeanFloralDiversity",
                 "rare.degree",
                 "(1|Site)")


## **********************************************************
## community model, check assumptions first before adding parasites
## **********************************************************

## run_plot_freq_model_diagnostics(remove_subset_formula(formula.flower.div),
##                                 this_data=spec.net[spec.net$Weights == 1,],
##                                 this_family="students", site.lat=site.or.lat)

## run_plot_freq_model_diagnostics(remove_subset_formula(formula.flower.abund),
##                                 this_data=spec.net[spec.net$Weights == 1,],
##                                 this_family="students",site.lat=site.or.lat)

## run_plot_freq_model_diagnostics(remove_subset_formula(formula.bee.abund),
##                                 this_data=spec.net[spec.net$Weights == 1,],
##                                 this_family="students", site.lat=site.or.lat)

## run_plot_freq_model_diagnostics(remove_subset_formula(formula.bombus.abund),
##                                 this_data=spec.net[spec.net$Weights == 1,],
##                                 this_family="students", site.lat=site.or.lat)

## run_plot_freq_model_diagnostics(remove_subset_formula(formula.HB.abund),
##                                 this_data=spec.net[spec.net$Weights == 1,],
##                                 this_family="students", site.lat=site.or.lat)

## run_plot_freq_model_diagnostics(remove_subset_formula(formula.bee.div),
##                                 this_data=spec.net[spec.net$Weights == 1,],
##                                 this_family="gaussian",
##                                 site.lat=site.or.lat)

## **********************************************************
## Parasite presence
## **********************************************************
## bombus
fit.bombus <- runCombinedParasiteModels(spec.bombus, species.group="bombus",
                                        parasite = c("CrithidiaPresence", "ApicystisSpp"),
                                        xvars=xvars.multi.bombus,
                                        ncores,
                                        iter = 10^4,
                                        chains = 1,
                                        thin=1,
                                        init=0, data2= list(phylo_matrix=phylo_matrix),
                                        SEM = TRUE, neg.binomial =  FALSE,
                                        site.lat=site.or.lat)


## melissodes
## there are not enough positives to get this model to converge
table(spec.melissodes$CrithidiaPresence[spec.melissodes$WeightsPar == 1])
table(spec.melissodes$ApicystisSpp[spec.melissodes$WeightsPar ==1 ])

fit.melissodes <- runCombinedParasiteModels(spec.melissodes, species.group="melissodes",
                                            parasite = c("CrithidiaPresence", "ApicystisSpp"),
                                            xvars=xvars.single.species,
                                            ncores,
                                            iter = 2*10^4,
                                            chains = 1,
                                            thin=1,
                                            init=0,
                                            SEM = TRUE,
                                            neg.binomial = TRUE,
                                            site.lat=site.or.lat)

## Honey bees
## there are not enough positives to get this model to converge
table(spec.apis$CrithidiaPresence[spec.melissodes$WeightsPar == 1])
table(spec.apis$ApicystisSpp[spec.melissodes$WeightsPar ==1 ])

fit.apis <- runCombinedParasiteModels(spec.apis, species.group="apis",
                                      parasite = c("CrithidiaPresence", "ApicystisSpp"),
                                      xvars=xvars.apis,
                                      ncores,
                                      iter = 4*10^4,
                                      chains = 1,
                                      thin=1,
                                      init=0,
                                      SEM = TRUE,
                                      neg.binomial = TRUE,
                                      site.lat=site.or.lat)



## heat maps

getParComm <- function(parasite, spec){
    parasite <- aggregate(list(Parasite=spec[, parasite]),
                          list(GenusSpecies=spec$GenusSpecies,
                               Site=spec$Site),
                          function(x) sum(x) / length(x))

    parasite.comm <- samp2site.spp(parasite$Site,
                                   parasite$GenusSpecies,
                                   parasite$Parasite)
    return(parasite.comm)
}


## screened
rownames(spec.net) <- NULL
spec.screened <- spec.net[spec.net$WeightsPar == 1,]
sp.n <- c(table(spec.screened$GenusSpecies))
sp.n <- sp.n[sp.n >=4]
spec.screened <- spec.screened[spec.screened$GenusSpecies %in% names(sp.n),]

spec.screened <- spec.screened[spec.screened$WeightsSp ==1,]

parasite.cols <- c( "SpCrithidiaPresence",
                   paste0("Sp", parasites))

spec.screened <- spec.screened[, c("GenusSpecies", "Genus", "Site",
                                   "SampleRound", "Year",
                                   "SpScreened",
                                   parasite.cols)]

sum.genus.screened <- spec.screened  %>%
    group_by(Site, Genus) %>%
    summarise(
        SpCrithidiaPresence= sum(SpCrithidiaPresence),
        SpApicystisSpp = sum(SpApicystisSpp),
        SpNosemaBombi= sum(SpNosemaBombi),
        SpNosemaCeranae = sum(SpNosemaCeranae),
        SpAscosphaeraSpp= sum(SpAscosphaeraSpp),
        SpCrithidiaExpoeki = sum(SpCrithidiaExpoeki),
        SpCrithidiaBombi = sum(SpCrithidiaBombi),
        SpCrithidiaSpp = sum(SpCrithidiaSpp),
        SpScreened = sum(SpScreened)
    )


sum.screened <- spec.screened  %>%
    group_by(Site, GenusSpecies) %>%
    summarise(
        SpCrithidiaPresence= sum(SpCrithidiaPresence),
        SpApicystisSpp = sum(SpApicystisSpp),
        SpNosemaBombi= sum(SpNosemaBombi),
        SpNosemaCeranae = sum(SpNosemaCeranae),
        SpAscosphaeraSpp= sum(SpAscosphaeraSpp),
        SpCrithidiaExpoeki = sum(SpCrithidiaExpoeki),
        SpCrithidiaBombi = sum(SpCrithidiaBombi),
        SpCrithidiaSpp = sum(SpCrithidiaSpp),
        SpScreened = sum(SpScreened)
        )

sum.screened[, parasite.cols ] <- sum.screened[, parasite.cols
                                               ]/sum.screened$SpScreened
sum.genus.screened[, parasite.cols ] <-
    sum.genus.screened[, parasite.cols ]/sum.genus.screened$SpScreened

sum.screened$SpScreened <- NULL
sum.genus.screened$SpScreened <- NULL



### need to update samp to site species
source("src/misc.R")
getParComm <- function(par.sum.dat, bee.col, par.name){
    parasite.comm <- samp2site.spp(par.sum.dat$Site,
                                   par.sum.dat[, bee.col],
                                   par.sum.dat[, par.name])
    return(parasite.comm)
}


getParComm(sum.screened, "GenusSpecies", "SpApicystisSpp")

## heat maps of # of infected individuals

plotParasiteMap <- function(){
    colfunc <- colorRampPalette(c("white", "red"))
    par(oma=c(8,4,3,2), mar=c(1,2,2,1),
        mgp=c(1.5,0.5,0))
    heatmap.2(parasite.comms[[parasite]],
              trace="none",
              col=colfunc,
              breaks=seq(0, 1, 0.1))
    mtext(parasite, 3, line=2.5)
}

parasite.comms <- lapply(parasites, getParComm)
names(parasite.comms) <- parasites



for(parasite in parasites){
    pdf.f(plotParasiteMap,
          file=sprintf("figures/heatmaps/%s.pdf", parasite),
          width=8, height=7)
}

