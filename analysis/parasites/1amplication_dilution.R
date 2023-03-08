setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')
## setwd('~/Dropbox (University of Oregon)/skyislands')

setwd("analysis/parasites")
rm(list=ls())
source("src/misc.R")
source("src/init.R")
source("src/writeResultsTable.R")
source("src/makeMultiLevelData.R")

ncores <- 2

## **********************************************************
## formula for site effects on the bee community
## **********************************************************
## order data by stand for ease of viewing
spec <- spec[order(spec$Site),]

spec$Lat <- log(spec$Lat)
spec$Area <- log(spec$Area)

## create a dumby varaible "Weight" to deal with the data sets being at
## different levels to get around the issue of having to pass in one
## data set into brms
spec$YearSR <- paste(spec$Year, spec$SampleRound, sep=";")

## will need to modify when we have multiple years
spec <- makeDataMultiLevel(spec, "Site", "YearSR")

## create a dumby varaible "WeightPar" for the parasite data. The
## original intention was to keep stan from dropping data for
## site-level models, but weight is 0 for parasite models.

spec$WeightsPar <- 1
spec$WeightsPar[is.na(spec$ParasitePresence) | spec$Apidae != 1] <- 0

## stan drops all NA data, so can set ParasitePresence to 0 with WeightsPar
## to keep it in the models
spec$ParasitePresence[is.na(spec$ParasitePresence | spec$Apidae != 1)] <- 0

## all of the variables that are explanatory variables and thus need
## to be centered
vars <- c("MeanFloralAbundance",
          "MeanFloralDiversity",
          "PollDiversity",
          "Lat", "SRDoy",
          "MeanITD",
          "r.degree")

spec.all <- spec
spec.all[, vars] <- apply(spec.all[, vars], 2, standardize)

## define all the formulas for the different parts of the models
## **********************************************************
## Model 1.1: formula for forest effects on floral community
## **********************************************************
## flower diversity
formula.flower.div <- formula(MeanFloralDiversity | weights(Weights) ~
                                  Lat +
                                      SRDoy + I(SRDoy^2) +
                                      (1|Site)
                              )
## flower abund
formula.flower.abund <- formula(MeanFloralAbundance | weights(Weights) ~
                                    SRDoy + I(SRDoy^2) +
                                        (1|Site)
                                )

## **********************************************************
## Model 1.2: formula for forest effects on bee community
## **********************************************************

## bee diversity
formula.bee.div <- formula(PollDiversity | weights(Weights)~
                               MeanFloralDiversity +
                                   Lat +  SRDoy +
                                   (1|Site)
                           )

## bombus abund
formula.bombus.abund <- formula(BombusAbundance | weights(Weights)~
                                    MeanFloralAbundance +
                                        MeanFloralDiversity+
                                        SRDoy + I(SRDoy^2) +
                                        Lat +
                                        (1|Site)
                                )
## HB abund
formula.HB.abund <- formula(HBAbundance | weights(Weights)~
                                MeanFloralAbundance +
                                    MeanFloralDiversity+
                                    SRDoy + I(SRDoy^2) +
                                    Lat +
                                    (1|Site)
                            )
## bee abund
formula.bee.abund <- formula(NonBombusHBAbundance | weights(Weights)~
                                 MeanFloralAbundance +
                                     MeanFloralDiversity+
                                     SRDoy + I(SRDoy^2) +
                                     Lat +
                                     (1|Site)
                             )

## **********************************************************
## Model 1.3: formula for bee community effects on parasitism
## **********************************************************
formula.parasite <- formula(ParasitePresence | weights(WeightsPar) ~
                                NonBombusHBAbundance +
                                    HBAbundance +
                                    BombusAbundance + Lat +
                                    PollDiversity + MeanITD + r.degree +
                                    (1|Site)
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
spec.bombus <- spec.all
spec.bombus$WeightsPar[spec.bombus$Genus != "Bombus"] <- 0

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
                                NonBombusHBAbundance +
                                    HBAbundance +
                                    BombusAbundance + Lat +
                                    PollDiversity +
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
spec.apis <- spec.all
spec.apis$WeightsPar[spec.apis$Genus != "Apis"] <- 0
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

## **********************************************************
## Model 1c parasite specific models
## **********************************************************
## Bombus - CrithidiaBombi
## **********************************************************
formula.CrithidiaBombi <- formula(CrithidiaBombi | weights(WeightsPar) ~
                                      NonBombusHBAbundance +
                                          HBAbundance +
                                          BombusAbundance + Lat +
                                          PollDiversity +
                                          (1|Site)
                                  )
bf.CrithidiaBombi <- bf(formula.CrithidiaBombi, family="bernoulli")
bform.CrithidiaBombi <- bf.fabund + bf.fdiv +
    bf.babund + bf.bombusabund + bf.HBabund +
    bf.bdiv + bf.CrithidiaBombi +
    set_rescor(FALSE)
fit.bombus.CrithidiaBombi <- brm(bform.CrithidiaBombi, spec.bombus,
                                 cores=ncores,
                                 iter = 10^4,
                                 chains = 1,
                                 thin=1,
                                 init=0,
                                 control = list(adapt_delta = 0.99))
write.ms.table(fit.bombus.CrithidiaBombi, "parasitism_bombus_CrithidiaBombi")
save(fit.bombus.CrithidiaBombi, spec.bombus,
     file="saved/parasiteFitBombus_CrithidiaBombi.Rdata")

## **********************************************************
## Bombus - CrithidiaExpoeki
## **********************************************************
formula.CrithidiaExpoeki <- formula(CrithidiaExpoeki | weights(WeightsPar) ~
                                      NonBombusHBAbundance +
                                          HBAbundance +
                                          BombusAbundance + Lat +
                                          PollDiversity +
                                          (1|Site)
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
## Bombus - NosemaCeranae
## **********************************************************
formula.NosemaCeranae <- formula(NosemaCeranae | weights(WeightsPar) ~
                                     NonBombusHBAbundance +
                                         HBAbundance +
                                         BombusAbundance + Lat +
                                         PollDiversity +
                                         (1|Site)
                                 )
bf.NosemaCeranae <- bf(formula.NosemaCeranae, family="bernoulli")
bform.NosemaCeranae <- bf.fabund + bf.fdiv +
    bf.babund + bf.bombusabund + bf.HBabund +
    bf.bdiv + bf.NosemaCeranae +
    set_rescor(FALSE)
fit.bombus.NosemaCeranae <- brm(bform.NosemaCeranae, spec.bombus,
                                cores=ncores,
                                iter = 10^4,
                                chains = 1,
                                thin=1,
                                init=0,
                                control = list(adapt_delta = 0.99))
write.ms.table(fit.bombus.NosemaCeranae,
               "parasitism_bombus_NosemaCeranae")
save(fit.bombus.NosemaCeranae, spec.bombus,
     file="saved/parasiteFitBombus_NosemaCeranae.Rdata")

## **********************************************************
## Model 1c - remove r.degree and body size from the models
## **********************************************************
## **********************************************************
## Apis - NosemaCeranae
## **********************************************************

formula.NosemaCeranae <- formula(NosemaCeranae | weights(WeightsPar) ~
                                     NonBombusHBAbundance +
                                         HBAbundance +
                                         BombusAbundance + Lat +
                                         PollDiversity +
                                         (1|Site)
                                 )
bf.NosemaCeranae <- bf(formula.parasite, family="bernoulli")
bform.NosemaCeranae <- bf.fabund + bf.fdiv +
    bf.babund + bf.bombusabund + bf.HBabund +
    bf.bdiv + bf.NosemaCeranae +
    set_rescor(FALSE)

fit.NosemaCeranae <- brm(bform.NosemaCeranae, spec.apis,
                cores=ncores,
                iter = 10^4,
                chains = 1,
                thin=1,
                init=0,
                control = list(adapt_delta = 0.99))

write.ms.table(fit.NosemaCeranae, "NosemaCeranae_apis")

save(fit.apis, spec,
     file="saved/parasiteFitApis_NosemaCeranae.Rdata")

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
