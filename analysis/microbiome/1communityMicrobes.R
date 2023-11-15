## setwd('C:/Users/na_ma/Dropbox (University of Oregon)/Rotation/skyIslands')
## setwd('~/Dropbox (University of Oregon)/skyIslands')
## setwd('~/Dropbox (University of Oregon)/skyIslands')
setwd("/Volumes/bombus/rhayes/Dropbox (University of Oregon)/skyIslands")

setwd("analysis/microbiome/")

rm(list=ls())
library(picante)
library(bayesplot)
library(brms)
source("src/misc_microbe.R")
source("src/misc.R")
source('src/makeMultiLevelData.R')
source("src/writeResultsTable.R")
source("src/runParasiteModels.R")
source("src/standardize_weights.R")
source("src/plant_poll_models_microbe.R")

load("../../data/spec_RBCL_16s.Rdata")

## all of the variables that are explanatory variables and thus need
## to be centered
vars_yearsr <- c("MeanFloralAbundance",
                 "MeanFloralDiversity",
                 "Net_BeeDiversity",
                 "Lat", "SRDoy"  
)
vars_sp <- c("MeanITD",
             "rare.degree")

variables.to.log <- "rare.degree"


source("src/init.R")
#source('src/init_microbe.R')


ncores <- 1






## all of the variables that are explanatory variables and thus need
## to be centered
vars_yearsr <- c("MeanFloralAbundance",
                 "MeanFloralDiversity",
                 "Net_BeeDiversity",
                 "Lat", "SRDoy"  
)
vars_sp <- c("MeanITD",
             "rare.degree")

variables.to.log <- "rare.degree"

## uses only net specimens, and drops syrphids
source("src/init_microbe.R")


#genus_pd_fit <- function(spec.net, this_genus, num_iter){
# 
# microbes <- colnames(spec.net)[grepl("16s:", colnames(spec.net))] 
# 
# screened.microbes <- apply(spec.net, 1, function(x) all(is.na(x[microbes])))
# 
# spec.microbes <- spec.net[!screened.microbes, ]
# 
# #genus.microbes <- spec.microbes[spec.microbes$Genus == this_genus, ]
# 
# 
# ## QUESTION: should include root = TRUE? if false gives warning 3x
# ## warning: Rooted tree and include.root=TRUE argument required to calculate PD of single-species communities. Single species community assigned PD value of NA.
# PD <- apply(spec.microbes[,microbes], 1, function(x){
#   this.bee <- x[x > 0]
#   this.tree <- prune.sample(t(this.bee), tree.16s)
#   pd(t(this.bee), this.tree, include.root = TRUE)
# })
# 
# PD <- do.call(rbind, PD)
# 
# spec.microbes <- cbind(spec.microbes, PD)
# 
# spec.net <- merge(spec.net, spec.microbes, all.x=TRUE)
# 
load("../../data/spec_RBCL_16s_PD.Rdata")

## define all the formulas for the different parts of the models
source("src/plant_poll_models_microbe.R")
spec.net <- spec.net.merge

## QUESTION: should there be NAs?
## check ids
unique(spec.net$GenusSpecies[spec.net$Apidae == 1 &
                               is.na(spec.net$MeanITD)])

# having trouble getting code to generate PD to run on lab computer, saved it on RH machine and loading it here

## **********************************************************
## Parasite models set up
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
                     iter = 10^4,
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
## Model 1.3: formula for bee community effects on parasitism
## **********************************************************

get_genus_PD_formula <- function(input.df, genus){
  genus.microbes <- input.df[input.df$Genus == genus, ]
  
  
  genus.PD <- apply(genus.microbes[,microbes], 1, function(x){
    this.bee <- x[x > 0]
    this.tree <- prune.sample(t(this.bee), tree.16s)
    pd(t(this.bee), this.tree, include.root = FALSE)
  })
  
  genus.PD <- do.call(rbind, genus.PD)
  
  genus.microbes <- cbind(genus.microbes, genus.PD)
  
  input.df <- merge(input.df, genus.microbes, all.x=TRUE)

  formula.microbes <- formula(genus.PD ~
                                PollAbundance*FloralDiversity +
                                PollDiversity +
                                FloralAbundance +
                                (1|Site)
  )
}

#Apis
apis.formula <- get_genus_PD_formula(spec, 'Apis')
#Bombus
bombus.formula <- get_genus_PD_formula(spec, 'Bombus')
#Megachile
megachile.formula <- get_genus_PD_formula(spec, 'Megachile')
#Anthophora
anthophora.formula <- get_genus_PD_formula(spec, 'Anthophora')



## **********************************************************
## Community models
## **********************************************************
## convert to brms format
bf.fabund <- bf(formula.flower.abund)
bf.fdiv <- bf(formula.flower.div)
bf.babund <- bf(formula.bee.abund)
bf.bdiv <- bf(formula.bee.div)
bf.apis <- bf(apis.formula)
bf.bombus <- bf(bombus.formula)
bf.megachile <- bf(megachile.formula)
bf.anthophora <- bf(anthophora.formula)

## **********************************************************
## Model 1 community effects on gut microbe phylo distance
## **********************************************************



#Apis model

apis_bform <- bf.fabund + bf.fdiv + bf.babund + bf.bdiv + bf.apis +
  set_rescor(FALSE)

## run model
apis_fit <- brm(apis_bform, spec,
           cores=ncores,
           iter = 10^3,
           chains = 2,
           thin=1,
           init=0,
           control = list(adapt_delta = 0.99))

write.ms.table(apis_fit, "apis")

#Bombus model
bombus_bform <- bf.fabund + bf.fdiv + bf.babund + bf.bdiv + bf.bombus +
  set_rescor(FALSE)

## run model
bombus_fit <- brm(bombus_bform, spec,
                cores=ncores,
                iter = 10^3,
                chains = 2,
                thin=1,
                init=0,
                control = list(adapt_delta = 0.99))


write.ms.table(bombus_fit, "bombus")

#Megachile

megachile_bform <- bf.fabund + bf.fdiv + bf.babund + bf.bdiv + bf.megachile +
  set_rescor(FALSE)

## run model
megachile_fit <- brm(megachile_bform, spec,
                cores=ncores,
                iter = 10^3,
                chains = 2,
                thin=1,
                init=0,
                control = list(adapt_delta = 0.99))

write.ms.table(megachile_fit, "megachile")

#Anthophora
anthophora_bform <- bf.fabund + bf.fdiv + bf.babund + bf.bdiv + bf.anthophora +
  set_rescor(FALSE)

## run model
anthophora_fit <- brm(anthophora_bform, spec,
                cores=ncores,
                iter = 10^3,
                chains = 2,
                thin=1,
                init=0,
                control = list(adapt_delta = 0.99))

write.ms.table(anthophora_fit, "anthophora")

#write.ms.table(fit, "microbes")


#save(fit, spec,
#     file="saved/microbesFitMod.Rdata")
## dignostic figures
#mcmc_trace(fit)
#ggsave("figures_diagnostics_microbes.pdf",
#       height=11, width=8.5)









