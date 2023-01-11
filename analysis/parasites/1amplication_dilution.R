setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')
## setwd('~/Dropbox (University of Oregon)/skyislands')

setwd("analysis/parasites")
rm(list=ls())
source("src/misc.R")
source("src/init.R")
source("src/writeResultsTable.R")
source("src/makeMultiLevelData.R")
ncores <- 10

## **********************************************************
## formula for site effects on the bee community
## **********************************************************
## order data by stand for ease of viewing
spec <- spec[order(spec$Site),]

## drop 2022 for now because not enough species IDed
spec <- spec[spec$Year != 2022,]
spec$Lat <- log(spec$Lat)
spec$Area <- log(spec$Area)

## all of the variables that are explanatory variables and thus need
## to be centered
vars <- c("MeanFloralAbundance",
          "MeanFloralDiversity",
          ## "PollAbundance",
          ## "HBAbundance",
          ## "BombusAbundance",
          ## "NonBombusHBAbundance",
          "PollDiversity",
          "Lat", "SRDoy",
          "Area",
          "MeanITD",
          "r.degree")

##  center all of the x variables across the datasets
spec[, vars] <- apply(spec[, vars], 2, standardize)

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
## to keep it in the models, this is commented out because we don't
## want to use the entire dataset in this publication

spec$ParasitePresence[is.na(spec$ParasitePresence | spec$Apidae != 1)] <- 0

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
## ## bee abund
## formula.bee.abund <- formula(PollAbundance | weights(Weights)~
##                                  MeanFloralAbundance +
##                                      SRDoy + I(SRDoy^2) +
##                                      (1|Site)
##                              )

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
                                    BombusAbundance +
                                    PollDiversity + MeanITD + r.degree +
                                    (1|Site)
                            )

## **********************************************************
## Community models
## **********************************************************
## convert to brms format
bf.fabund <- bf(formula.flower.abund)
bf.fdiv <- bf(formula.flower.div)
bf.babund <- bf(formula.bee.abund, family="Poisson")
bf.bombusabund <- bf(formula.bombus.abund, family="Poisson")
bf.HBabund <- bf(formula.HB.abund, family="Poisson")
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
fit <- brm(bform, spec,
           cores=ncores,
           iter = 10^4,
           chains = 1,
           thin=1,
           init=0,
           control = list(adapt_delta = 0.99))
write.ms.table(fit, "parasitism")
save(fit, spec,
     file="saved/parasiteFitMod.Rdata")
## dignostic figures
plot.res(fit, "parasite")


## social only
fit.social <- brm(bform, spec[spec$Genus == "Bombus" |
                              spec$Genus == "Apis",],
                  cores=ncores,
                  iter = 10^4,
                  chains = 1,
                  thin=1,
                  init=0,
                  control = list(adapt_delta = 0.99))

write.ms.table(fit.social, "parasitism_social")

save(fit.social, spec,
     file="saved/parasiteSocialFitMod.Rdata")



## Model checks
mod.floral <- lmer(MeanFloralDiversity ~   Lat + Area +   (1|Site),
                   data=spec[spec$Weights == 1,])
vif(mod.floral)

## including elevation is very colinear

mod.floral.abund <- lmer(MeanFloralAbundance ~ Area + Lat +  (1|Site),
                         data=spec[spec$Weights == 1,])
vif(mod.floral.abund)
