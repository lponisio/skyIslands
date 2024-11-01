## **********************************************************
## Load libraries
## **********************************************************

rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path
setwd(local.path)
setwd("skyIslands/analysis/microbiome/")


run.diagnostics = FALSE
make.plots = FALSE
run.bombus = TRUE
run.apis = FALSE
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


## **********************************************************
## Standardize, center, and transform data
## **********************************************************

## all of the variables that are explanatory variables and thus need
## to be centered

variables.to.log <- c("rare.degree",
                      "BeeAbundance"
)

vars_yearsr <- c("MeanFloralAbundance",
                 "MeanFloralDiversity",
                 "SRDoy",
                 "BeeAbundance",
                 "BeeDiversity",
                 "VisitedFloralDiversity"
)

vars_yearsrsp <- "rare.degree"
vars_sp <- "MeanITD.x"
vars_site <- "Lat"

## **********************************************************
## Source files
## **********************************************************

source("src/misc_microbe.R")
source("src/misc.R")
source('src/makeMultiLevelData.R')
source("src/standardize_weights_parasites.R") #TODO fix title
source("src/init_microbe.R")
source("src/writeResultsTable.R")
source("src/runPlotFreqModelDiagnostics.R")


ncores <- 1


## **********************************************************
## Flower diversity
## **********************************************************

flower.div.vars <- c("Lat",
                     "(1|Site)")

flower.div.x <- paste(flower.div.vars, collapse="+")
flower.div.y <- "MeanFloralDiversity | subset(Weights)"
formula.flower.div <- as.formula(paste(flower.div.y, "~",flower.div.x))

## **********************************************************
## Bee abundance
## **********************************************************

tot.bee.abund.vars <- c("(1|Site)")

tot.bee.abund.x <- paste(tot.bee.abund.vars, collapse="+")
tot.bee.abund.y <- "BeeAbundance | subset(Weights)"
formula.tot.bee.abund <- as.formula(paste(tot.bee.abund.y, "~",tot.bee.abund.x))

## **********************************************************
## Bee diversity
## **********************************************************

tot.bee.div.vars <- c("MeanFloralDiversity",
                      "Lat",
                      "(1|Site)")

tot.bee.div.x <- paste(tot.bee.div.vars, collapse="+")
tot.bee.div.y <- "BeeDiversity | subset(Weights)"
formula.tot.bee.div <- as.formula(paste(tot.bee.div.y, "~",tot.bee.div.x))

## **********************************************************
## convert formulas to brms format
## **********************************************************
bf.fdiv <- bf(formula.flower.div)
bf.tot.babund <- bf(formula.tot.bee.abund)
bf.bdiv <- bf(formula.bee.div)
bf.tot.bdiv <- bf(formula.tot.bee.div)

## **********************************************************
## add weights for each genus and microbe type
## **********************************************************

## TODO add to init if works
spec.net$WeightsObligateBombus = spec.net$WeightsObligateMicrobe*spec.net$BombusWeights

## TODO add to init if works
spec.net$WeightsTransientBombus = spec.net$WeightsTransientMicrobe*spec.net$BombusWeights

## TODO add to init if works

# use pd obligate with skew normal
spec.net$WeightsObligateApis = spec.net$WeightsObligateMicrobe*spec.net$ApisWeights

#use pd transient with skew normal?
spec.net$WeightsTransientApis = spec.net$WeightsTransientMicrobe*spec.net$ApisWeights

# use pd obligate lof with skew normal?
spec.net$WeightsObligateMelissodes = spec.net$WeightsObligateMicrobe*spec.net$MelissodesWeights
# use pd transient log with gaussian or student t?
spec.net$WeightsTransientMelissodes = spec.net$WeightsTransientMicrobe*spec.net$MelissodesWeights


## **********************************************************
## Bombus microbe PD model
## **********************************************************

## Obligate model
ob.microbe.bombus.vars <- c("BeeAbundance",
                         "BeeDiversity",
                         "Lat", 
                         "MeanFloralDiversity",
                         "MeanITD.x",
                         "rare.degree",
                         "(1|Site)",
                         "(1|gr(GenusSpecies, cov = phylo_matrix))") 


ob.microbe.bombus.x <- paste(ob.microbe.bombus.vars, collapse="+")
ob.microbe.bombus.y <- "PD.obligate.log | subset(WeightsObligateBombus)"
formula.ob.microbe.bombus <- as.formula(paste(ob.microbe.bombus.y, "~",
                                           ob.microbe.bombus.x))

bf.ob.microbe.bombus.skew <- bf(formula.ob.microbe.bombus, family=skew_normal())

## Facultative model
non.ob.microbe.bombus.vars <- c("BeeAbundance",
                            "BeeDiversity",
                            "Lat", 
                            "MeanFloralDiversity",
                            "MeanITD.x",
                            "rare.degree",
                            "(1|Site)",
                            "(1|gr(GenusSpecies, cov = phylo_matrix))") 

non.ob.microbe.bombus.x <- paste(non.ob.microbe.bombus.vars, collapse="+")
non.ob.microbe.bombus.y <- "PD.transient.log | subset(WeightsTransientBombus)"
formula.non.ob.microbe.bombus <- as.formula(paste(non.ob.microbe.bombus.y, "~",
                                              non.ob.microbe.bombus.x))

bf.non.ob.microbe.bombus.student <- bf(formula.non.ob.microbe.bombus, family=student())


## combined model
bform.bombus <- bf.fdiv +
  bf.tot.bdiv +
  bf.tot.babund +
  bf.ob.microbe.bombus.skew +
  bf.non.ob.microbe.bombus.student +
  set_rescor(FALSE)

## TODO check if need to change all NAs to 0
spec.net[is.na(spec.net)] <- 0

if(run.bombus){
  fit.microbe.bombus <- brm(bform.bombus, spec.net,
                          cores=ncores,
                          iter = 10000,
                          chains =1,
                          thin=1,
                          init=0,
                          open_progress = FALSE,
                          control = list(adapt_delta = 0.99,
                                         max_treedepth = 15),
                          save_pars = save_pars(all = TRUE),
                          data2 = list(phylo_matrix=phylo_matrix))
  
  write.ms.table(fit.microbe.bombus, "bombus_microbe_2")
  r2loo.bombus <- loo_R2(fit.microbe.bombus)
  r2.bombus <- rstantools::bayes_R2(fit.microbe.bombus)
  save(fit.microbe.bombus, spec.net, r2.bombus, r2loo.bombus,
       file="saved/fullMicrobeBombusFit_2.Rdata")
}


## **********************************************************
## Apis microbe PD model
## **********************************************************
## obligate PD model
ob.microbe.apis.vars <- c("BeeAbundance",
                          "BeeDiversity",
                          "Lat", 
                          "MeanFloralDiversity",
                          "rare.degree",
                            "(1|Site)") 


ob.microbe.apis.x <- paste(ob.microbe.apis.vars, collapse="+")
ob.microbe.apis.y <- "PD.obligate | subset(WeightsObligateApis)"
formula.ob.microbe.apis <- as.formula(paste(ob.microbe.apis.y, "~",
                                              ob.microbe.apis.x))

## set up brms formulas for different family types
bf.ob.microbe.apis.skew <- bf(formula.ob.microbe.apis, family=skew_normal())
bf.ob.microbe.apis.student <- bf(formula.ob.microbe.apis, family=student())
bf.ob.microbe.apis.gaussian <- bf(formula.ob.microbe.apis)
bf.ob.microbe.apis.mix.sg <- bf(formula.ob.microbe.apis, family=mixture(student, gaussian))
bf.ob.microbe.apis.mix.gg <- bf(formula.ob.microbe.apis, family=mixture(gaussian, gaussian))
bf.ob.microbe.apis.mix.ss <- bf(formula.ob.microbe.apis, family=mixture(student, student))



## facultative microbe pd model
non.ob.microbe.apis.vars <- c("BeeAbundance",
                              "BeeDiversity",
                              "Lat",
                              "MeanFloralDiversity",
                              "rare.degree",
                              "(1|Site)") 


non.ob.microbe.apis.x <- paste(non.ob.microbe.apis.vars, collapse="+")
non.ob.microbe.apis.y <- "PD.transient | subset(WeightsTransientApis)"
formula.non.ob.microbe.apis <- as.formula(paste(non.ob.microbe.apis.y, "~",
                                                  non.ob.microbe.apis.x))


## set up brms formulas for different family types
bf.non.ob.microbe.apis.skew <- bf(formula.non.ob.microbe.apis, family=skew_normal())
bf.non.ob.microbe.apis.student <- bf(formula.non.ob.microbe.apis, family=student())
bf.non.ob.microbe.apis.gaussian <- bf(formula.non.ob.microbe.apis)



## combined model

#combine forms
bform.apis <- bf.fdiv +
   bf.tot.bdiv +
   bf.tot.babund +
  bf.ob.microbe.apis.student +
  bf.non.ob.microbe.apis.student +
  set_rescor(FALSE)

if(run.apis){
fit.microbe.apis <- brm(bform.apis , spec.net,
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
save(fit.microbe.apis, spec.net, r2.apis, r2loo.apis,
     file="saved/fullMicrobeApisFit.Rdata")
}

## **********************************************************
## Melissodes microbe PD model
## **********************************************************
## obligate PD model
ob.microbe.melissodes.vars <- c("BeeAbundance",
                          "BeeDiversity",
                          "Lat", 
                          "MeanFloralDiversity",
                          "rare.degree",
                          "(1|Site)") 


ob.microbe.melissodes.x <- paste(ob.microbe.melissodes.vars, collapse="+")
ob.microbe.melissodes.y <- "PD.obligate.log | subset(WeightsObligateMelissodes)"
formula.ob.microbe.melissodes <- as.formula(paste(ob.microbe.melissodes.y, "~",
                                            ob.microbe.melissodes.x))

## set up model families in brms format
bf.ob.microbe.melissodes.skew <- bf(formula.ob.microbe.melissodes, family=skew_normal())
bf.ob.microbe.melissodes.student <- bf(formula.ob.microbe.melissodes, family=student())
bf.ob.microbe.melissodes.gaussian <- bf(formula.ob.microbe.melissodes)


## facultative PD model
non.ob.microbe.melissodes.vars <- c("BeeAbundance",
                              "BeeDiversity",
                              "Lat", 
                              "MeanFloralDiversity",
                              "rare.degree",
                              "(1|Site)") 


non.ob.microbe.melissodes.x <- paste(non.ob.microbe.melissodes.vars, collapse="+")
non.ob.microbe.melissodes.y <- "PD.transient.log | subset(WeightsTransientMelissodes)"
formula.non.ob.microbe.melissodes <- as.formula(paste(non.ob.microbe.melissodes.y, "~",
                                                non.ob.microbe.melissodes.x))


## set up model families in brms format
bf.non.ob.microbe.melissodes.student <- bf(formula.non.ob.microbe.melissodes, family=student())
bf.non.ob.microbe.melissodes.skew <- bf(formula.non.ob.microbe.melissodes, family=skew_normal())
bf.non.ob.microbe.melissodes.gaussian <- bf(formula.non.ob.microbe.melissodes)


##combine forms
bform.melissodes <- bf.fdiv +
  bf.tot.bdiv +
  bf.tot.babund +
  bf.ob.microbe.melissodes.skew +
  bf.non.ob.microbe.melissodes.skew +
  set_rescor(FALSE)

if(run.melissodes){
fit.microbe.melissodes <- brm(bform.melissodes , spec.net,
                           cores=ncores,
                           iter = 10000,
                           chains =1,
                           thin=1,
                           init=0,
                           open_progress = FALSE,
                           control = list(adapt_delta = 0.9999,
                                          max_treedepth = 15
                                          ),
                           save_pars = save_pars(all = TRUE))

write.ms.table(fit.microbe.melissodes, "melissodes_microbe_2")
r2loo.melissodes <- loo_R2(fit.microbe.melissodes)
r2.melissodes <- rstantools::bayes_R2(fit.microbe.melissodes)
save(fit.microbe.melissodes, spec.net, r2.melissodes, r2loo.melissodes,
      file="saved/fullMicrobeMelissodesFit_2.Rdata")
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