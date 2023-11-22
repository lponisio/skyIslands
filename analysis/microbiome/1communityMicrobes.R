## setwd('C:/Users/na_ma/Dropbox (University of Oregon)/Rotation/skyIslands')
## setwd('~/Dropbox (University of Oregon)/skyIslands')
setwd('~/Dropbox (University of Oregon)/skyIslands')
#setwd("/Volumes/bombus/rhayes/Dropbox (University of Oregon)/skyIslands")

setwd("analysis/microbiome/")

rm(list=ls())

run.diagnostics = TRUE

library(picante)
library(bayesplot)
library(pscl)
library(brms)
library(performance)
library(lme4)
library(glmmADMB)
library(R2admb)
# install.packages("glmmADMB", 
#                  repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                          getOption("repos")),
#                  type="source")

## all of the variables that are explanatory variables and thus need
## to be centered
vars_yearsr <- c("MeanFloralAbundance",
                 "MeanFloralDiversity",
                 "Net_BeeDiversity",
                 "Lat", "SRDoy", "Net_BombusAbundance", "Net_HBAbundance"  
)
vars_sp <- c("MeanITD",
             "rare.degree")

variables.to.log <- c("rare.degree",
                      "MeanFloralAbundance")

source("src/misc_microbe.R")
source("src/misc.R")
source('src/makeMultiLevelData.R')
source("src/standardize_weights.R")
source("src/init_microbe.R")
source("src/writeResultsTable.R")
source("src/runParasiteModels.R")
source("src/runPlotFreqModelDiagnostics.R")


ncores <- 1


## QUESTION: should include root = TRUE? if false gives warning 3x
## warning: Rooted tree and include.root=TRUE argument required to calculate PD of single-species communities. Single species community assigned PD value of NA.
PD <- apply(spec.microbes[,microbes], 1, function(x){
  this.bee <- x[x > 0]
  this.tree <- prune.sample(t(this.bee), tree.16s)
  #browser()
  picante::pd(t(this.bee), this.tree, include.root = TRUE)
})

PD <- do.call(rbind, PD)

spec.microbes <- cbind(spec.microbes, PD)

spec.net <- merge(spec.net, spec.microbes, all.x=TRUE)



## QUESTION: should there be NAs? not sure what check ids means
## check ids
unique(spec.net$GenusSpecies[spec.net$Apidae == 1 &
                               is.na(spec.net$MeanITD)])

#dropping genus species with NA -- 
#spec.net <- spec.net[!is.na(spec.net$GenusSpecies),]



## **********************************************************
## Flower abundance
## **********************************************************

## flower abundance variables 
flower.abund.vars <- c("Year",
                       "SRDoy",
                       "I(SRDoy^2)",
                       "Lat",
                       "(1|Site)")

flower.abund.x <- paste(flower.abund.vars, collapse="+")
flower.abund.y <- "MeanFloralAbundance | weights(Weights)"
formula.flower.abund <- as.formula(paste(flower.abund.y, "~",flower.abund.x))
                                         

#VegAbund check
freq.formula.flower.abund <- as.formula(paste("MeanFloralAbundance", "~", flower.abund.x ))

#for this_data, use spec.net[spec.net$Weights==1,] to incorporate weights into frequentist models
if(run.diagnostics){
  freq.model.flower.abund <- run_plot_freq_model_diagnostics(
    freq.formula.flower.abund,
    spec.net[spec.net$Weights==1,],
    this_family = 'gaussian')
  
  ggsave(freq.model.flower.abund,
         file="figures/SI_VegAbundModelDiagnostics.pdf",
         height=8, width=11)
}

## **********************************************************
## Flower diversity
## **********************************************************


## flower abundance variables 
flower.div.vars <- c("Year",
                       "SRDoy",
                       "I(SRDoy^2)",
                       "Lat",
                       "(1|Site)")

flower.div.x <- paste(flower.div.vars, collapse="+")
flower.div.y <- "MeanFloralDiversity | weights(Weights)"
formula.flower.div <- as.formula(paste(flower.div.y, "~",flower.div.x))


#Vegdiv check
freq.formula.flower.div <- as.formula(paste("MeanFloralDiversity", "~", flower.div.x ))

#for this_data, use spec.net[spec.net$Weights==1,] to incorporate weights into frequentist models
if(run.diagnostics){
  freq.model.flower.div <- run_plot_freq_model_diagnostics(
    freq.formula.flower.div,
    spec.net[spec.net$Weights==1,],
    this_family = 'gaussian')
  
  ggsave(freq.model.flower.div,
         file="figures/SI_VegDivModelDiagnostics.pdf",
         height=8, width=11)
}

## **********************************************************
## Bee abundance
## **********************************************************

## bombus abundance variables 
bombus.abund.vars <- c("MeanFloralAbundance",
                       "Year",
                       "SRDoy",
                       "I(SRDoy^2)",
                       "Lat",
                       "(1|Site)")

bombus.abund.x <- paste(bombus.abund.vars, collapse="+")
bombus.abund.y <- "Net_BombusAbundance | weights(Weights)"
formula.bombus.abund <- as.formula(paste(bombus.abund.y, "~",bombus.abund.x))


#Bombus abundance check
freq.formula.bombus.abund <- as.formula(paste("Net_BombusAbundance", "~", bombus.abund.x ))

#for this_data, use spec.net[spec.net$Weights==1,] to incorporate weights into frequentist models
if(run.diagnostics){
  freq.model.bombus.abund <- run_plot_freq_model_diagnostics(
    freq.formula.bombus.abund,
    spec.net[spec.net$Weights==1,],
    this_family = "zero_inflated_negbinomial")
  
  ggsave(freq.model.bombus.abund,
         file="figures/SI_BombusAbundModelDiagnostics.pdf",
         height=8, width=11)
}


## HB abund

## honeybee abundance variables 
hb.abund.vars <- c("MeanFloralAbundance",
                       "Year",
                       "SRDoy",
                       "I(SRDoy^2)",
                       "Lat",
                       "(1|Site)")

hb.abund.x <- paste(hb.abund.vars, collapse="+")
hb.abund.y <- "Net_HBAbundance | weights(Weights)"
formula.hb.abund <- as.formula(paste(hb.abund.y, "~",hb.abund.x))


# HB abund check
freq.formula.hb.abund <- as.formula(paste("Net_HBAbundance", "~", hb.abund.x ))

#for this_data, use spec.net[spec.net$Weights==1,] to incorporate weights into frequentist models
if(run.diagnostics){
  freq.model.hb.abund <- run_plot_freq_model_diagnostics(
    freq.formula.hb.abund,
    spec.net[spec.net$Weights==1,],
    this_family = "negbinomial")
  
  ggsave(freq.model.hb.abund,
         file="figures/SI_HoneyBeeAbundModelDiagnostics.pdf",
         height=8, width=11)
}


## bee abund

bee.abund.vars <- c("MeanFloralAbundance",
                   "Year",
                   "SRDoy",
                   "I(SRDoy^2)",
                   "Lat",
                   "(1|Site)")

bee.abund.x <- paste(bee.abund.vars, collapse="+")
bee.abund.y <- "Net_NonBombusHBAbundance | weights(Weights)"
formula.bee.abund <- as.formula(paste(bee.abund.y, "~",bee.abund.x))


# bee abund check
freq.formula.bee.abund <- as.formula(paste("Net_NonBombusHBAbundance", "~", bee.abund.x ))

##for this_data, use spec.net[spec.net$Weights==1,] to incorporate weights into frequentist models
if(run.diagnostics){
  freq.model.bee.abund <- run_plot_freq_model_diagnostics(
    freq.formula.bee.abund,
    spec.net[spec.net$Weights==1,],
    this_family = "zero_inflated_negbinomial")
  
  ggsave(freq.model.bee.abund,
         file="figures/SI_BeeAbundModelDiagnostics.pdf",
         height=8, width=11)
}


## Warning message:
##   In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
##                  Model failed to converge with max|grad| = 0.00307909 (tol = 0.002, component 1)


## **********************************************************
## Bee diversity
## **********************************************************


bee.div.vars <- c("MeanFloralDiversity",
                    "Year",
                    "SRDoy",
                    "I(SRDoy^2)",
                    "Lat",
                    "(1|Site)")

bee.div.x <- paste(bee.div.vars, collapse="+")
bee.div.y <- "Net_BeeDiversity | weights(Weights)"
formula.bee.div <- as.formula(paste(bee.div.y, "~",bee.div.x))


# bee div check
freq.formula.bee.div <- as.formula(paste("Net_BeeDiversity", "~", bee.div.x ))

##for this_data, use spec.net[spec.net$Weights==1,] to incorporate weights into frequentist models
if(run.diagnostics){
  freq.model.bee.div <- run_plot_freq_model_diagnostics(
    freq.formula.bee.div,
    spec.net[spec.net$Weights==1,],
    this_family = "gaussian")
  
  ggsave(freq.model.bee.div,
         file="figures/SI_BeeDiversityModelDiagnostics.pdf",
         height=8, width=11)
}


## **********************************************************
## Microbe models set up
## **********************************************************
## Multi species models
xvars.multi.species <-  c("Net_NonBombusHBAbundance",
                          "Net_HBAbundance",
                          "Net_BombusAbundance",
                          "Net_BeeDiversity",
                          "rare.degree", "MeanITD",
                          "(1|Site)", "(1|GenusSpecies)")
## single species models
xvars.single.species <-  c("Net_NonBombusHBAbundance",
                           "Net_HBAbundance",
                           "Net_BombusAbundance",
                           "Net_BeeDiversity",
                           "rare.degree",
                           "(1|Site)")

## **********************************************************
## community model
## **********************************************************

bform.community <- bf.fabund + bf.fdiv +
  bf.babund +
  bf.bombusabund + bf.HBabund +
  bf.bdiv  +
  set_rescor(FALSE)

fit.community <- brm(bform.community, spec.net,
                     cores=ncores,
                     iter = 100,
                     chains = 1,
                     thin=1,
                     init=0,
                     control = list(adapt_delta = 0.99),
                     save_pars = save_pars(all = TRUE))
write.ms.table(fit.community,
               sprintf("microbe_%s_%s",
                       species.group="all", parasite="none"))
r2loo <- loo_R2(fit.community)
r2 <- rstantools::bayes_R2(fit.community)
save(fit.community, spec.net, r2,
     file="saved/communityFit.Rdata")


## **********************************************************
## Model 1 community effects on gut microbe phylo distance
## **********************************************************


microbe.vars <- xvars.multi.species

microbe.x <- paste(microbe.vars, collapse="+")
microbe.y <- "PD | weights(WeightsPar)"
formula.microbe <- as.formula(paste(microbe.y, "~",
                                     microbe.x))


bf.microbe <- bf(formula.microbe)

#combine forms
bform <- bf.fabund + bf.fdiv +
  bf.babund + bf.bombusabund + bf.HBabund +
  bf.bdiv  +    
  bf.microbe +
  set_rescor(FALSE)

## run model
fit.microbe <- brm(bform , spec.microbes,
                  cores=ncores,
                  iter = 10000,
                  chains =1,
                  thin=1,
                  init=0,
                  open_progress = FALSE,
                  control = list(adapt_delta = 0.99),
                  save_pars = save_pars(all = TRUE))

write.ms.table(fit.microbe,
               sprintf("full_microbe_%s_%s",
                       species.group="all", parasite="none"))
r2loo <- loo_R2(fit.microbe)
r2 <- rstantools::bayes_R2(fit.microbe)
save(fit.microbe, spec.net, r2, r2loo,
     file="saved/fullMicrobeFit.Rdata")
