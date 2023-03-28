setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')
## setwd('~/Dropbox (University of Oregon)/skyislands')

setwd("analysis/parasites")
rm(list=ls())
source("src/misc.R")
source("src/writeResultsTable.R")
source("src/makeMultiLevelData.R")

## uses only net specimens, and drops syrphids
source("src/init.R")
ncores <- 1

## **********************************************************
## formula for site effects on the bee community
## **********************************************************
## order data by stand for ease of viewing
spec.net <- spec.net[order(spec.net$Site),]

spec.net$Lat <- log(spec.net$Lat)
spec.net$Area <- log(spec.net$Area)

## create a dumby varaible "Weight" to deal with the data sets being at
## different levels to get around the issue of having to pass in one
## data set into brms
spec.net$YearSR <- paste(spec.net$Year, spec.net$SampleRound, sep=";")

## will need to modify when we have multiple years
spec.net <- makeDataMultiLevel(spec.net, "Site", "YearSR")

## create a dumby varaible "WeightPar" for the parasite data. The
## original intention was to keep stan from dropping data for
## site-level models, but weight is 0 for parasite models.

spec.net$WeightsPar <- 1
spec.net$WeightsPar[is.na(spec.net$ParasitePresence) | spec.net$Apidae != 1] <- 0

## stan drops all NA data, so can set ParasitePresence to 0 with WeightsPar
## to keep it in the models
spec.net$ParasitePresence[is.na(spec.net$ParasitePresence | spec.net$Apidae != 1)] <- 0

unique(spec.net$GenusSpecies[spec.net$Apidae == 1 & is.na(spec.net$MeanITD)])

## all of the variables that are explanatory variables and thus need
## to be centered
vars <- c("MeanFloralAbundance",
          "MeanFloralDiversity",
          "PollDiversity",
          "Lat", "SRDoy",
           "MeanITD",
          ## "r.degree", ## across all networks
          "rare.degree"  ## site-year level
          )

spec.all <- spec.net
spec.all[, vars] <- apply(spec.all[, vars], 2, standardize)

## bombus only data
spec.bombus <- spec.all
spec.bombus$WeightsPar[spec.bombus$Genus != "Bombus"] <- 0

## apis only data
spec.apis <- spec.all
spec.apis$WeightsPar[spec.apis$Genus != "Apis"] <- 0

## melissodes only data
spec.melissodes <- spec.all
spec.melissodes$WeightsPar[spec.melissodes$Genus != "Melissodes"] <- 0

spec.apidae <- spec.all
spec.apidae$WeightsPar[spec.melissodes$Family != "Apidae"] <- 0

## screened table for model planning
screened <- spec.net[spec.net$Apidae == 1,]
screened.tab <- table(screened$GenusSpecies, screened$Apidae, screened$Site)
screened.tab

## not enough replication of species in Halictidae
spec.all$lWeightsPar[spec.all$Family =="Halictidae"] <- 0
## not enough replication of species in Colletidae
spec.all$lWeightsPar[spec.all$Family =="Colletidae"] <- 0

## define all the formulas for the different parts of the models
## **********************************************************
## Model 1.1: formula for forest effects on floral community
## **********************************************************
## flower diversity
formula.flower.div <- formula(MeanFloralDiversity | weights(Weights) ~
                                  Lat +  SRDoy + I(SRDoy^2) + Year +
                                      (1|Site)
                              )
## flower abund
formula.flower.abund <- formula(MeanFloralAbundance | weights(Weights) ~
                                    SRDoy + I(SRDoy^2) + Year +
                                        (1|Site)
                                )

## **********************************************************
## Model 1.2: formula for forest effects on bee community
## **********************************************************
## bee diversity
formula.bee.div <- formula(Net_BeeDiversity | weights(Weights)~
                               MeanFloralDiversity +
                                   Lat +  SRDoy + Year +
                                   (1|Site)
                           )
## bombus abund
formula.bombus.abund <- formula(Net_BombusAbundance | weights(Weights)~
                                    MeanFloralAbundance +
                                        MeanFloralDiversity+
                                        SRDoy + I(SRDoy^2) +
                                        Lat + Year+
                                        (1|Site)
                                )
## HB abund
formula.HB.abund <- formula(Net_HBAbundance | weights(Weights)~
                                MeanFloralAbundance +
                                    MeanFloralDiversity+
                                    SRDoy + I(SRDoy^2) +
                                    Lat + Year+
                                    (1|Site)
                            )
## bee abund
formula.bee.abund <- formula(Net_NonBombusHBAbundance | weights(Weights)~
                                 MeanFloralAbundance +
                                     MeanFloralDiversity+
                                     SRDoy + I(SRDoy^2) +
                                     Lat + Year+
                                     (1|Site)
                             )

## **********************************************************
## Model 1.3: formula for bee community effects on parasitism
## **********************************************************
formula.parasite <- formula(ParasitePresence | weights(WeightsPar) ~
                              Net_NonBombusHBAbundance +
                                Net_HBAbundance +
                                Net_BombusAbundance +
                                Net_BeeDiversity + MeanITD + rare.degree +
                                (1|Site) + (1|GenusSpecies) + (1|Genus)
                            )

## **********************************************************
## Community models
## **********************************************************
## convert to brms format
bf.fabund <- bf(formula.flower.abund)
bf.fdiv <- bf(formula.flower.div)
bf.babund <- bf(formula.bee.abund, family="negbinomial")
bf.bombusabund <- bf(formula.bombus.abund, family="negbinomial")
bf.HBabund <- bf(formula.HB.abund, family="negbinomial")
bf.bdiv <- bf(formula.bee.div)
bf.par <- bf(formula.parasite, family="bernoulli")

## **********************************************************
## Model 1 community effects on bee parasitism
## **********************************************************
## full model
bform <- bf.fabund + bf.fdiv +
    bf.babund + bf.bombusabund + bf.HBabund +
    bf.bdiv + bf.par +
    set_rescor(FALSE)

## run model
fit <- brm(bform, spec.all,
           cores=ncores,
           iter = 10^4,
           chains = 1,
           thin=1,
           init=0,
           control = list(adapt_delta = 0.99))
write.ms.table(fit, "parasitism")
save(fit, spec.all,
     file="saved/parasiteFitMod.Rdata")
## dignostic figures
plot.res(fit, "parasite")

## **********************************************************
## Model 1b community effects on bee parasitism
## **********************************************************
## bombus only
## **********************************************************
fit.bombus <- brm(bform, spec.bombus,
                  cores=ncores,
                  iter = 10^4,
                  chains = 1,
                  thin=1,
                  init=0,
                  control = list(adapt_delta = 0.99))

write.ms.table(fit.bombus, "parasitism_bombus")

save(fit.bombus, spec.bombus,
     file="saved/parasiteFitBombusMod.Rdata")

## **********************************************************
## Model 1b ## remove r.degree and body size from the models
## **********************************************************
## Apis and Melissodes
## **********************************************************
formula.parasite <- formula(ParasitePresence | weights(WeightsPar) ~
                              Net_NonBombusHBAbundance +
                                Net_HBAbundance +
                                Net_BombusAbundance +
                                Net_BeeDiversity +
                                rare.degree
                                (1|Site)
                            )
bf.par <- bf(formula.parasite, family="bernoulli")
bform <- bf.fabund + bf.fdiv +
  bf.babund + bf.bombusabund + bf.HBabund +
  bf.bdiv + bf.par +
  set_rescor(FALSE)

## **********************************************************
## Apis only
## **********************************************************
fit.apis <- brm(bform, spec.apis,
                cores=ncores,
                iter = 10^4,
                chains = 1,
                thin=1,
                init=0,
                control = list(adapt_delta = 0.99))
write.ms.table(fit.apis, "parasitism_apis")
save(fit.apis, spec,
     file="saved/parasiteFitApisMod.Rdata")

## **********************************************************
## Melissodes only
## **********************************************************
fit.melissodes <- brm(bform, spec.melissodes,
                      cores=ncores,
                      iter = 10^4,
                      chains = 1,
                      thin=1,
                      init=0,
                      control = list(adapt_delta = 0.99))
write.ms.table(fit.melissodes, "parasitism_melissodes")
save(fit.melissodes, spec.melissodes,
     file="saved/parasiteFitMelissodesMod.Rdata")

## **********************************************************
## Model 1c parasite specific models
## **********************************************************
## Bombus - CrithidiaBombi
## **********************************************************

xvars.multi.species <-  c("Net_NonBombusHBAbundance",
                          "Net_HBAbundance",
                          "Net_BombusAbundance",
                          "Net_BeeDiversity",
                          "rare.degree", "MeanITD",
                          "(1|Site)", "(1|GenusSpecies)")

runParasiteModels <- function(spec.data,
                              species.group, parasite,
                              xvars){

    formula.parasite  <- as.formula(paste(
        paste(parasite, "| weights(WeightsPar)"),
        paste(xvars,
              collapse=" + "),
        sep=" ~ "))


    bf.parasite <- bf(formula.parasite, family="bernoulli")
    bform.parasite <- bf.fabund + bf.fdiv +
        bf.babund + bf.bombusabund + bf.HBabund +
        bf.bdiv + bf.parasite +
        set_rescor(FALSE)

    fit.parasite <- brm(bform.parasite, spec.data,
                        cores=ncores,
                        iter = 10^4,
                        chains = 1,
                        thin=1,
                        init=0,
                        control = list(adapt_delta = 0.99))
    write.ms.table(fit.parasite,
                   sprintf("parasitism_%s_%s",
                           species.group, parasite))
    save(fit.parasite, spec.data,
         file=sprintf("saved/parasiteFit_%s_%s.Rdata",
                      species.group, parasite))
    return(fit.parasite)
}

fit.bombus.cbombi <- runParasiteModelMultiSpecies(spec.bombus, "bombus",
                                                  "CrithidiaBombi",
                                                  xvars.multi.species)

fit.bombus.cexpoeki <- runParasiteModelMultiSpecies(spec.bombus, "bombus",
                                                  "CrithidiaExpoeki",
                                                  xvars.multi.species)

fit.bombus.cspp <- runParasiteModelMultiSpecies(spec.bombus, "bombus",
                                                  "CrithidiaPresence",
                                                  xvars.multi.species)

## **********************************************************
## Bombus - CrithidiaExpoeki
## **********************************************************
formula.CrithidiaExpoeki <- formula(CrithidiaExpoeki | weights(WeightsPar) ~
                                      NonBombusHBAbundance +
                                          HBAbundance +
                                          BombusAbundance +
                                          Net_BeeDiversity +
                                         rare.degree + MeanITD +
                                          (1|Site) + (1|GenusSpecies)
                                  )
bf.CrithidiaExpoeki <- bf(formula.CrithidiaExpoeki, family="bernoulli")
bform.CrithidiaExpoeki <- bf.fabund + bf.fdiv +
    bf.babund + bf.bombusabund + bf.HBabund +
    bf.bdiv + bf.CrithidiaExpoeki +
    set_rescor(FALSE)
fit.bombus.CrithidiaExpoeki <- brm(bform.CrithidiaExpoeki, spec.bombus,
                                 cores=ncores,
                                 iter = 10^4,
                                 chains = 1,
                                 thin=1,
                                 init=0,
                                 control = list(adapt_delta = 0.99))
write.ms.table(fit.bombus.CrithidiaExpoeki,
               "parasitism_bombus_CrithidiaExpoeki")
save(fit.bombus.CrithidiaExpoeki, spec.bombus,
     file="saved/parasiteFitBombus_CrithidiaExpoeki.Rdata")

## **********************************************************
## Bombus - Any crithidia
## **********************************************************
formula.CrithidiaPresence <- formula(CrithidiaPresence | weights(WeightsPar) ~
                                      NonBombusHBAbundance +
                                          HBAbundance +
                                          BombusAbundance +
                                          Net_BeeDiversity +
                                          rare.degree + MeanITD +
                                          (1|Site) + (1|GenusSpecies)
                                  )
bf.CrithidiaPresence <- bf(formula.CrithidiaPresence, family="bernoulli")
bform.CrithidiaPresence <- bf.fabund + bf.fdiv +
    bf.babund + bf.bombusabund + bf.HBabund +
    bf.bdiv + bf.CrithidiaPresence +
    set_rescor(FALSE)
fit.bombus.CrithidiaPresence <- brm(bform.CrithidiaPresence, spec.bombus,
                                 cores=ncores,
                                 iter = 10^4,
                                 chains = 1,
                                 thin=1,
                                 init=0,
                                 control = list(adapt_delta = 0.99))
write.ms.table(fit.bombus.CrithidiaPresence,
               "parasitism_bombus_CrithidiaPresence")
save(fit.bombus.CrithidiaPresence, spec.bombus,
     file="saved/parasiteFitBombus_CrithidiaPresence.Rdata")

## ## **********************************************************
## ## Bombus - NosemaCeranae
## ## **********************************************************
## formula.NosemaCeranae <- formula(NosemaCeranae | weights(WeightsPar) ~
##                                      NonBombusHBAbundance +
##                                          HBAbundance +
##                                          BombusAbundance + Lat +
##                                          Net_BeeDiversity +
##                                          (1|Site)
##                                  )
## bf.NosemaCeranae <- bf(formula.NosemaCeranae, family="bernoulli")
## bform.NosemaCeranae <- bf.fabund + bf.fdiv +
##     bf.babund + bf.bombusabund + bf.HBabund +
##     bf.bdiv + bf.NosemaCeranae +
##     set_rescor(FALSE)
## fit.bombus.NosemaCeranae <- brm(bform.NosemaCeranae, spec.bombus,
##                                 cores=ncores,
##                                 iter = 10^4,
##                                 chains = 1,
##                                 thin=1,
##                                 init=0,
##                                 control = list(adapt_delta = 0.99))
## write.ms.table(fit.bombus.NosemaCeranae,
##                "parasitism_bombus_NosemaCeranae")
## save(fit.bombus.NosemaCeranae, spec.bombus,
##      file="saved/parasiteFitBombus_NosemaCeranae.Rdata")

## **********************************************************
## Model 1c - remove r.degree and body size from the models
## **********************************************************

## ## *******************************************************
## ## Apis - NosemaCeranae
## ## *******************************************************
## formula.NosemaCeranae <- formula(NosemaCeranae | weights(WeightsPar) ~
##                                      NonBombusHBAbundance +
##                                          HBAbundance +
##                                          BombusAbundance + Lat +
##                                          Net_BeeDiversity +
##                                          (1|Site)
##                                  )
## bf.NosemaCeranae <- bf(formula.parasite, family="bernoulli")
## bform.NosemaCeranae <- bf.fabund + bf.fdiv +
##     bf.babund + bf.bombusabund + bf.HBabund +
##     bf.bdiv + bf.NosemaCeranae +
##     set_rescor(FALSE)

## fit.NosemaCeranae <- brm(bform.NosemaCeranae, spec.apis,
##                 cores=ncores,
##                 iter = 10^4,
##                 chains = 1,
##                 thin=1,
##                 init=0,
##                 control = list(adapt_delta = 0.99))

## write.ms.table(fit.NosemaCeranae, "NosemaCeranae_apis")

## save(fit.apis, spec,
##      file="saved/parasiteFitApis_NosemaCeranae.Rdata")

## **********************************************************
## Melissodes only
## **********************************************************
spec.melissodes <- spec.all
spec.melissodes$WeightsPar[spec.melissodes$Genus != "Melissodes"] <- 0
fit.melissodes <- brm(bform, spec.melissodes,
                      cores=ncores,
                      iter = 10^4,
                      chains = 1,
                      thin=1,
                      init=0,
                      control = list(adapt_delta = 0.99))

write.ms.table(fit.melissodes, "parasitism_melissodes")

save(fit.melissodes, spec.melissodes,
     file="saved/parasiteFitMelissodesMod.Rdata")
