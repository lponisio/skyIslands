rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path
setwd(local.path)
setwd("skyIslands/analysis/microbiome/")


run.diagnostics = FALSE
make.plots = FALSE
run.bombus = TRUE
run.apis = TRUE
run.melissodes = TRUE

library(picante)
library(plyr)
library(bayesplot)
library(pscl)
library(brms)
library(performance)
#library(lme4)
#library(glmmADMB)
library(R2admb)
library(shinystan)
# install.packages("glmmADMB", 
#                  repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                          getOption("repos")),
#                  type="source")

## all of the variables that are explanatory variables and thus need
## to be centered

##QUESTION: log or center first?? 
# 
variables.to.log <- c("rare.degree",
                      "MeanFloralAbundance", #keep logged
                      "BeeAbundance"
)

vars_yearsr <- c("MeanFloralAbundance",
                 "MeanFloralDiversity",
                 "SRDoy",
                 "BeeAbundance",
                 "BeeDiversity",
                 "VisitedFloralDiversity"
                 #"FloralDiversity"
)

vars_yearsrsp <- "rare.degree"
vars_sp <- "MeanITD"
vars_site <- "Lat"



source("src/misc_microbe.R")
source("src/misc.R")
source('src/makeMultiLevelData.R')
source("src/standardize_weights_parasites.R")
# incorporating Nicoles standardize functions
#source("src/standardize_weights.R")
source("src/init_microbe.R")
source("src/writeResultsTable.R")
source("src/runPlotFreqModelDiagnostics.R")


ncores <- 1


## **********************************************************
## Flower abundance
## **********************************************************

## flower abundance variables 
flower.abund.vars <- c("Year",
                       #"SRDoy",
                       #"I(SRDoy^2)",
                       "(1|Site)")

flower.abund.x <- paste(flower.abund.vars, collapse="+")
flower.abund.y <- "MeanFloralAbundance | subset(Weights)"
formula.flower.abund <- as.formula(paste(flower.abund.y, "~",flower.abund.x))




## **********************************************************
## Flower diversity
## **********************************************************


## flower abundance variables 
flower.div.vars <- c("Year",
                     #"SRDoy",
                     #"I(SRDoy^2)",
                     "Lat",
                     "(1|Site)")

flower.div.x <- paste(flower.div.vars, collapse="+")
flower.div.y <- "MeanFloralDiversity | subset(Weights)"
formula.flower.div <- as.formula(paste(flower.div.y, "~",flower.div.x))



## **********************************************************
## Bee abundance
## **********************************************************



## bee abund total
tot.bee.abund.vars <- c("MeanFloralAbundance",
                        "Year",
                        #"SRDoy",
                        #"I(SRDoy^2)",
                        "(1|Site)")

tot.bee.abund.x <- paste(tot.bee.abund.vars, collapse="+")
tot.bee.abund.y <- "BeeAbundance | subset(Weights)"
formula.tot.bee.abund <- as.formula(paste(tot.bee.abund.y, "~",tot.bee.abund.x))




#net bee abund
## bee abund total
net.bee.abund.vars <- c("MeanFloralAbundance",
                        "Year",
                        #"SRDoy",
                        #"I(SRDoy^2)",
                        "(1|Site)")

net.bee.abund.x <- paste(net.bee.abund.vars, collapse="+")
net.bee.abund.y <- "Net_BeeAbundance | subset(Weights)"
formula.net.bee.abund <- as.formula(paste(net.bee.abund.y, "~",net.bee.abund.x))





## **********************************************************
## Bee diversity
## **********************************************************


bee.div.vars <- c("MeanFloralDiversity",
                  "Year",
                  #"SRDoy",
                  #"I(SRDoy^2)",
                  "Lat",
                  "(1|Site)")

bee.div.x <- paste(bee.div.vars, collapse="+")
bee.div.y <- "Net_BeeDiversity | subset(Weights)"
formula.bee.div <- as.formula(paste(bee.div.y, "~",bee.div.x))

## bee div total
tot.bee.div.vars <- c("MeanFloralDiversity",
                      "Year",
                      #"SRDoy",
                      #"I(SRDoy^2)",
                      "Lat",
                      "(1|Site)")

tot.bee.div.x <- paste(tot.bee.div.vars, collapse="+")
tot.bee.div.y <- "BeeDiversity | subset(Weights)"
formula.tot.bee.div <- as.formula(paste(tot.bee.div.y, "~",tot.bee.div.x))



## **********************************************************
## convert formulas to brms forma
## **********************************************************
bf.fabund <- bf(formula.flower.abund)
bf.fdiv <- bf(formula.flower.div)
bf.tot.babund <- bf(formula.tot.bee.abund)
bf.bdiv <- bf(formula.bee.div)
bf.tot.bdiv <- bf(formula.tot.bee.div)

## **********************************************************
## community model
## **********************************************************
# 
bform.community <- bf.fdiv +
  #bf.fabund +    #removing from mods because no hypotheses regarding fabund
  bf.tot.babund +
  bf.tot.bdiv  +
  set_rescor(FALSE)

# fit.community <- brm(bform.community, spec.net,
#                      cores=ncores,
#                      iter = 10000,
#                      chains = 1,
#                      thin=1,
#                      init=0,
#                      control = list(adapt_delta = 0.99),
#                      save_pars = save_pars(all = TRUE))
# write.ms.table(fit.community,
#                sprintf("community_%s_%s",
#                        species.group="all", parasite="none"))
# r2loo <- loo_R2(fit.community)
# r2 <- rstantools::bayes_R2(fit.community)
# save(fit.community, spec.net, r2,
#      file="saved/communityFit.Rdata")


## **********************************************************
## Model 1 community effects on gut microbe phylo distance
## **********************************************************

## **********************************************************
## Microbe models set up
## **********************************************************
## Multi species models

#cant put in sociality bc all ones with weightsPar == 1 are eusocial or missing

##probs want to add some measure of diet.. is mean plant diversity enough?
## not sure if we have RBCL richness yet

#
# microbe.vars <-  c("BeeAbundance",
#                     "BeeDiversity", "Lat", #check this doesn't make VIF high
#                     "MeanFloralDiversity", "MeanITD", "Sociality", # if not at the genus level
#                     "(1|Site)", "rare.degree", "(1|gr(GenusSpecies, cov = phylo_matrix))") # add cov matrix for each genus
# #rare.degree
# # # split genus into separate PD for apis bombus megachile
# #
# #
# microbe.x <- paste(microbe.vars, collapse="+")
# microbe.y <- "PD | weights(LogWeightsAbund)"
# formula.microbe <- as.formula(paste(microbe.y, "~",
#                                      microbe.x))
# #
# #
# bf.microbe <- bf(formula.microbe, family='lognormal')
# 
# # #combine forms
# bform <- bf.fabund +
#    bf.fdiv +
#    bf.tot.babund +
#    bf.tot.bdiv  +
#    bf.microbe +
#    set_rescor(FALSE)
#
# ## run model
#  fit.microbe <- brm(bform , spec.net,
#                     cores=ncores,
#                     iter = 5000,
#                     chains =1,
#                     thin=1,
#                     init=0,
#                     open_progress = FALSE,
#                     control = list(adapt_delta = 0.99),
#                     save_pars = save_pars(all = TRUE),
#                     data2 = list(phylo_matrix=phylo_matrix))
# 
# write.ms.table(fit.microbe, "full_microbe")
# r2loo <- loo_R2(fit.microbe)
# r2 <- rstantools::bayes_R2(fit.microbe)
# save(fit.microbe, spec.net, r2, r2loo,
#       file="saved/fullMicrobeFit.Rdata")

## run bombus model
microbe.bombus.vars <- c("BeeAbundance",
                                    "BeeDiversity", "Lat", #check this doesn't make VIF high
                                    "MeanFloralDiversity", "MeanITD",  "rare.degree",
                                    "(1|Site)", "(1|gr(GenusSpecies, cov = phylo_matrix))") # add cov matrix for each genus
## NA check
#check_for_NA(microbe.bombus.vars)


microbe.bombus.x <- paste(microbe.bombus.vars, collapse="+")
microbe.bombus.y <- "PD | subset(LogWeightsAbund)"
formula.microbe.bombus <- as.formula(paste(microbe.bombus.y, "~",
                                      microbe.bombus.x))


bf.microbe.bombus <- bf(formula.microbe.bombus)

## obligate PD model
ob.microbe.bombus.vars <- c("BeeAbundance",
                         "BeeDiversity", "Lat", #check this doesn't make VIF high
                         "MeanFloralDiversity", "MeanITD",  "rare.degree",
                         "(1|Site)", "(1|gr(GenusSpecies, cov = phylo_matrix))") # add cov matrix for each genus
## NA check
#check_for_NA(ob.microbe.bombus.vars)


ob.microbe.bombus.x <- paste(ob.microbe.bombus.vars, collapse="+")
ob.microbe.bombus.y <- "PD.obligate.log | subset(LogWeightsObligateAbund)"
formula.ob.microbe.bombus <- as.formula(paste(ob.microbe.bombus.y, "~",
                                           ob.microbe.bombus.x))

bf.ob.microbe.bombus.skew <- bf(formula.ob.microbe.bombus, family=skew_normal())
bf.ob.microbe.bombus.student <- bf(formula.ob.microbe.bombus, family=student())
bf.ob.microbe.bombus.g <- bf(formula.ob.microbe.bombus)

## non ob PD model
non.ob.microbe.bombus.vars <- c("BeeAbundance",
                            "BeeDiversity", "Lat", #check this doesn't make VIF high
                            "MeanFloralDiversity", "MeanITD",  "rare.degree",
                            "(1|Site)", "(1|gr(GenusSpecies, cov = phylo_matrix))") # add cov matrix for each genus
## NA check
check_for_NA(non.ob.microbe.bombus.vars)

non.ob.microbe.bombus.x <- paste(non.ob.microbe.bombus.vars, collapse="+")
non.ob.microbe.bombus.y <- "PD.transient.log | subset(LogWeightsTransientAbund)"
formula.non.ob.microbe.bombus <- as.formula(paste(non.ob.microbe.bombus.y, "~",
                                              non.ob.microbe.bombus.x))



bf.non.ob.microbe.bombus.skew <- bf(formula.non.ob.microbe.bombus, family=skew_normal())
bf.non.ob.microbe.bombus.student <- bf(formula.non.ob.microbe.bombus, family=student())
bf.non.ob.microbe.bombus.g <- bf(formula.non.ob.microbe.bombus)

#combine forms
bform.bombus.skew <- bf.ob.microbe.bombus.skew +
    bf.non.ob.microbe.bombus.skew +
    set_rescor(FALSE)

#combine forms
bform.bombus.student <- bf.ob.microbe.bombus.student +
  bf.non.ob.microbe.bombus.student +
  set_rescor(FALSE)
#combine forms
bform.bombus.g <- bf.ob.microbe.bombus.g +
  bf.non.ob.microbe.bombus.g +
  set_rescor(FALSE)

## combined model
#combine forms
bform.bombus <- bf.tot.babund +
  #bf.fdiv +
  #bf.tot.bdiv  +
  bf.ob.microbe.bombus.skew +
  bf.non.ob.microbe.bombus.student +
  set_rescor(FALSE)

## currently running with jsut the subset i tested the different microbe distribution fams with
## high rhats (2.07) and max treedepth issues 4487

## rerun with the full data

## 8/21/24 running mods with spec.bombus and full mod, now trying just microbs and bee diversity

## 8/21/24 mods with microb and bee div were going to take weeks to run, probs misspecified
##  trying just microbes and b abund

## change all NAs to 0
spec.bombus[is.na(spec.bombus)] <- 0

if(run.bombus){
  fit.microbe.bombus <- brm(bform.bombus, spec.bombus,
                          cores=ncores,
                          iter = 10000,
                          chains =1,
                          thin=1,
                          init=0,
                          open_progress = FALSE,
                          control = list(adapt_delta = 0.99),
                          save_pars = save_pars(all = TRUE),
                          data2 = list(phylo_matrix=phylo_matrix))
  
  write.ms.table(fit.microbe.bombus, "bombus_microbe")
  r2loo.bombus <- loo_R2(fit.microbe.bombus)
  r2.bombus <- rstantools::bayes_R2(fit.microbe.bombus)
  save(fit.microbe.bombus, spec.bombus, r2.bombus, r2loo.bombus,
       file="saved/fullMicrobeBombusFit.Rdata")
}

## rerunning models with just PD layers and dataset filtered to exclude any 0s in PD
spec.bombus2 <- spec.bombus[spec.bombus$PD.obligate.log != 0,]
spec.bombus2 <- spec.bombus[spec.bombus$PD.transient.log != 0,]


if(run.bombus){
  ## skew mod
fit.microbe.bombus.skew <- brm(bform.bombus.skew , spec.bombus2,
                     cores=ncores,
                      iter = 10000,
                     chains = 1,
                     thin=1,
                     init=0,
                     open_progress = FALSE,
                     control = list(adapt_delta = 0.99,
                                    max_treedepth=12),
                     save_pars = save_pars(all = TRUE),
                     data2 = list(phylo_matrix=phylo_matrix))

write.ms.table(fit.microbe.bombus.skew, "bombus_microbe_skew")
#r2loo.bombus <- loo_R2(fit.microbe.bombus)
r2.bombus.skew <- rstantools::bayes_R2(fit.microbe.bombus.skew)
save(fit.microbe.bombus.skew, spec.bombus2, r2.bombus.skew, #r2loo.bombus,
       file="saved/fullMicrobeBombusFit_skew.Rdata")

## student mod
fit.microbe.bombus.student <- brm(bform.bombus.student , spec.bombus2,
                               cores=ncores,
                               iter = 10000,
                               chains = 1,
                               thin=1,
                               init=0,
                               open_progress = FALSE,
                               control = list(adapt_delta = 0.99,
                                              max_treedepth=12),
                               save_pars = save_pars(all = TRUE),
                               data2 = list(phylo_matrix=phylo_matrix))

write.ms.table(fit.microbe.bombus.student, "bombus_microbe_student")
#r2loo.bombus <- loo_R2(fit.microbe.bombus)
r2.bombus.student <- rstantools::bayes_R2(fit.microbe.bombus.student)
save(fit.microbe.bombus.student, spec.bombus2, r2.bombus.student, #r2loo.bombus,
     file="saved/fullMicrobeBombusFit_student.Rdata")

## gaussian mod
fit.microbe.bombus.g <- brm(bform.bombus.g , spec.bombus2,
                                  cores=ncores,
                                  iter = 10000,
                                  chains = 1,
                                  thin=1,
                                  init=0,
                                  open_progress = FALSE,
                                  control = list(adapt_delta = 0.99,
                                                 max_treedepth=12),
                                  save_pars = save_pars(all = TRUE),
                                  data2 = list(phylo_matrix=phylo_matrix))

write.ms.table(fit.microbe.bombus.g, "bombus_microbe_g")
#r2loo.bombus <- loo_R2(fit.microbe.bombus)
r2.bombus.g <- rstantools::bayes_R2(fit.microbe.bombus.g)
save(fit.microbe.bombus.g, spec.bombus2, r2.bombus.g, #r2loo.bombus,
     file="saved/fullMicrobeBombusFit_g.Rdata")

}

# model comparison
# Assuming your models are named model1, model2, model3
waic1 <- waic(fit.microbe.bombus.g, resp='PDobligatelog')
waic2 <- waic(fit.microbe.bombus.skew, resp='PDobligatelog')
waic3 <- waic(fit.microbe.bombus.student, resp='PDobligatelog')

# Compare the models
loo_compare(waic1, waic2, waic3)



update.bombus = FALSE
if(update.bombus == TRUE){
  fit.microbe.bombus <- update(fit.microbe.bombus,
                               iter=30000,
                               control=list(max_treedepth=15))
  write.ms.table(fit.microbe.bombus, "bombus_microbe")
  r2loo.bombus <- loo_R2(fit.microbe.bombus)
  r2.bombus <- rstantools::bayes_R2(fit.microbe.bombus)
  save(fit.microbe.bombus, spec.bombus, r2.bombus, r2loo.bombus,
       file="saved/fullMicrobeBombusFit.Rdata")
}
## run apis model
microbe.apis.vars <- c("BeeAbundance",
                       "BeeDiversity", "Lat", #check this doesn't make VIF high
                       "MeanFloralDiversity",# "MeanITD",
                       "(1|Site)", "rare.degree"#, "(1|gr(GenusSpecies, cov = phylo_matrix))") # add cov matrix for each genus
)


microbe.apis.x <- paste(microbe.apis.vars, collapse="+")
microbe.apis.y <- "PD | weights(LogWeightsAbund)"
formula.microbe.apis <- as.formula(paste(microbe.apis.y, "~",
                                         microbe.apis.x))


bf.microbe.apis <- bf(formula.microbe.apis)

## run apis model
microbe.apis.vars <- c("BeeAbundance",
                       "BeeDiversity", "Lat", #check this doesn't make VIF high
                       "MeanFloralDiversity",# "MeanITD",
                       "(1|Site)", "rare.degree"#, "(1|gr(GenusSpecies, cov = phylo_matrix))") # add cov matrix for each genus
)


ob.microbe.apis.x <- paste(microbe.apis.vars, collapse="+")
ob.microbe.apis.y <- "PD.obligate.log | weights(LogWeightsObligateAbund)"
formula.ob.microbe.apis <- as.formula(paste(ob.microbe.apis.y, "~",
                                         ob.microbe.apis.x))


bf.ob.microbe.apis <- bf(formula.ob.microbe.apis)

# non obligate

non.ob.microbe.apis.x <- paste(microbe.apis.vars, collapse="+")
non.ob.microbe.apis.y <- "PD.transient.log | weights(LogWeightsTransientAbund)"
formula.non.ob.microbe.apis <- as.formula(paste(non.ob.microbe.apis.y, "~",
                                            non.ob.microbe.apis.x))


bf.non.ob.microbe.apis <- bf(formula.non.ob.microbe.apis)

#combine forms
bform.apis <- bf.fdiv +
  bf.tot.babund +
  bf.tot.bdiv  +
  bf.ob.microbe.apis +
  bf.non.ob.microbe.apis +
  set_rescor(FALSE)

if(run.apis){
fit.microbe.apis <- brm(bform.apis , spec.apis,
                        cores=ncores,
                        iter = 10000,
                        chains =1,
                        thin=1,
                        init=0,
                        open_progress = FALSE,
                        control = list(adapt_delta = 0.99),
                        save_pars = save_pars(all = TRUE))

write.ms.table(fit.microbe.apis, "apis_microbe")
r2loo.apis <- loo_R2(fit.microbe.apis)
r2.apis <- rstantools::bayes_R2(fit.microbe.apis)
save(fit.microbe.apis, spec.apis, r2.apis, r2loo.apis,
     file="saved/fullMicrobeApisFit.Rdata")
}
# ## run melissodes model
microbe.melissodes.vars <- c("BeeAbundance",
                         "BeeDiversity", "Lat", #check this doesn't make VIF high
                         "MeanFloralDiversity", #"MeanITD",
                         "(1|Site)", "rare.degree") #, "(1|gr(GenusSpecies, cov = phylo_matrix))") # add cov matrix for each genus



microbe.melissodes.x <- paste(microbe.melissodes.vars, collapse="+")
microbe.melissodes.y <- "PD | weights(LogWeightsAbund)"
formula.microbe.melissodes <- as.formula(paste(microbe.melissodes.y, "~",
                                           microbe.melissodes.x))


bf.microbe.melissodes <- bf(formula.microbe.melissodes)

## obligate


ob.microbe.melissodes.x <- paste(microbe.melissodes.vars, collapse="+")
ob.microbe.melissodes.y <- "PD.obligate.log | weights(LogWeightsObligateAbund)"
formula.ob.microbe.melissodes <- as.formula(paste(ob.microbe.melissodes.y, "~",
                                               ob.microbe.melissodes.x))


bf.ob.microbe.melissodes <- bf(formula.ob.microbe.melissodes)

## non obligate


non.ob.microbe.melissodes.x <- paste(microbe.melissodes.vars, collapse="+")
non.ob.microbe.melissodes.y <- "PD.transient.log | weights(LogWeightsTransientAbund)"
formula.non.ob.microbe.melissodes <- as.formula(paste(non.ob.microbe.melissodes.y, "~",
                                                  non.ob.microbe.melissodes.x))


bf.non.ob.microbe.melissodes <- bf(formula.non.ob.microbe.melissodes)

#combine forms
bform.melissodes <- bf.fdiv +
  bf.tot.babund +
  bf.tot.bdiv  +
  bf.ob.microbe.melissodes +
  bf.non.ob.microbe.melissodes +
  set_rescor(FALSE)

if(run.melissodes){
fit.microbe.melissodes <- brm(bform.melissodes , spec.melissodes,
                           cores=ncores,
                           iter = 10000,
                           chains =1,
                           thin=1,
                           init=0,
                           open_progress = FALSE,
                           control = list(adapt_delta = 0.99),
                           save_pars = save_pars(all = TRUE))

write.ms.table(fit.microbe.melissodes, "melissodes_microbe")
r2loo.melissodes <- loo_R2(fit.microbe.melissodes)
r2.melissodes <- rstantools::bayes_R2(fit.microbe.melissodes)
save(fit.microbe.melissodes, spec.melissodes, r2.melissodes, r2loo.melissodes,
      file="saved/fullMicrobeMelissodesFit.Rdata")
}


## **********************************************************
## Diagnostics
## **********************************************************
#VegAbund check
if(run.diagnostics){
  # freq.formula.flower.abund <- as.formula(paste("MeanFloralAbundance", "~", flower.abund.x ))
  # 
  # #for this_data, use spec.net[spec.net$Weights==1,] to incorporate weights into frequentist models
  # if(run.diagnostics){
  #   freq.model.flower.abund <- run_plot_freq_model_diagnostics(
  #     freq.formula.flower.abund,
  #     spec.net[spec.net$Weights==1,],
  #     this_family = 'gaussian')
  #   
  #   ggsave(freq.model.flower.abund,
  #          file="figures/diagnostics/SI_VegAbundModelDiagnostics.pdf",
  #          height=8, width=11)
  # }
  # 
  # #Vegdiv check
  # freq.formula.flower.div <- as.formula(paste("MeanFloralDiversity", "~", flower.div.x ))
  # 
  # #for this_data, use spec.net[spec.net$Weights==1,] to incorporate weights into frequentist models
  # if(run.diagnostics){
  #   freq.model.flower.div <- run_plot_freq_model_diagnostics(
  #     freq.formula.flower.div,
  #     spec.net[spec.net$Weights==1,],
  #     this_family = 'gaussian')
  #   
  #   ggsave(freq.model.flower.div,
  #          file="figures/diagnostics/SI_VegDivModelDiagnostics.pdf",
  #          height=8, width=11)
  # }
  # 
  # # bee abund check
  # freq.formula.tot.bee.abund <- as.formula(paste("BeeAbundance", "~", tot.bee.abund.x ))
  # 
  # ##for this_data, use spec.net[spec.net$Weights==1,] to incorporate weights into frequentist models
  # if(run.diagnostics){
  #   freq.model.tot.bee.abund <- run_plot_freq_model_diagnostics(
  #     freq.formula.tot.bee.abund,
  #     spec.net[spec.net$Weights==1,],
  #     this_family = "gaussian")
  #   
  #   ggsave(freq.model.tot.bee.abund,
  #          file="figures/diagnostics/SI_TotalBeeAbundModelDiagnostics.pdf",
  #          height=8, width=11)
  # }
  # 
  # # # bee abund check
  # # freq.formula.net.bee.abund <- as.formula(paste("Net_BeeAbundance", "~", net.bee.abund.x ))
  # # 
  # # ##for this_data, use spec.net[spec.net$Weights==1,] to incorporate weights into frequentist models
  # # if(run.diagnostics){
  # #   freq.model.net.bee.abund <- run_plot_freq_model_diagnostics(
  # #     freq.formula.net.bee.abund,
  # #     spec.net[spec.net$Weights==1,],
  # #     this_family = "gaussian")
  # #   
  # #   ggsave(freq.model.net.bee.abund,
  # #          file="figures/diagnostics/SI_NetBeeAbundModelDiagnostics.pdf",
  # #          height=8, width=11)
  # # }
  # 
  # # # bee div check
  # # freq.formula.bee.div <- as.formula(paste("Net_BeeDiversity", "~", bee.div.x ))
  # # 
  # # ##for this_data, use spec.net[spec.net$Weights==1,] to incorporate weights into frequentist models
  # # if(run.diagnostics){
  # #   freq.model.bee.div <- run_plot_freq_model_diagnostics(
  # #     freq.formula.bee.div,
  # #     spec.net[spec.net$Weights==1,],
  # #     this_family = "gaussian")
  # #   
  # #   ggsave(freq.model.bee.div,
  # #          file="figures/diagnostics/SI_BeeDiversityModelDiagnostics.pdf",
  # #          height=8, width=11)
  # # }
  # 
  # # total bee div check
  # freq.formula.tot.bee.div <- as.formula(paste("BeeDiversity", "~", tot.bee.div.x ))
  # 
  # ##for this_data, use spec.net[spec.net$Weights==1,] to incorporate weights into frequentist models
  # if(run.diagnostics){
  #   freq.model.tot.bee.div <- run_plot_freq_model_diagnostics(
  #     freq.formula.tot.bee.div,
  #     spec.net[spec.net$Weights==1,],
  #     this_family = "gaussian")
  #   
  #   ggsave(freq.model.tot.bee.div,
  #          file="figures/diagnostics/SI_TotalBeeDivModelDiagnostics.pdf",
  #          height=8, width=11)
  # }
  
  # microbe check
  freq.formula.microbe <- as.formula(paste("PD", "~", microbe.x ))
  
  ##for this_data, use spec.net[spec.net$Weights==1,] to incorporate weights into frequentist models
  if(run.diagnostics){
    freq.model.microbe <- run_plot_freq_model_diagnostics(
      freq.formula.microbe,
      spec.net[spec.net$WeightsPar==1,],
      this_family = "gaussian",
      launch.shiny = FALSE,
      examine.pairs = FALSE)
    
    ggsave(freq.model.microbe,
           file="figures/diagnostics/SI_MicrobeModelDiagnostics.pdf",
           height=8, width=11)
  }
}