## setwd('C:/Users/na_ma/Dropbox (University of Oregon)/Rotation/skyIslands')
## setwd('~/Dropbox (University of Oregon)/skyIslands')
setwd('~/Dropbox (University of Oregon)/skyIslands')
#setwd("/Volumes/bombus/rhayes/Dropbox (University of Oregon)/skyIslands")

setwd("analysis/microbiome/")

rm(list=ls())

run.diagnostics = FALSE

library(picante)
library(plyr)
library(bayesplot)
library(pscl)
library(brms)
library(performance)
library(lme4)
library(glmmADMB)
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
                 "Lat", "SRDoy",
                 "BeeAbundance",
                 "BeeDiversity",
                 "VisitedFloralDiversity"
                 #"FloralDiversity"
)

vars_sp <- c("MeanITD",
             "rare.degree")



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
# PD <- apply(spec.microbes[,microbes], 1, function(x){
#   this.bee <- x[x > 0]
#   this.tree <- prune.sample(t(this.bee), tree.16s)
#   #browser()
#   picante::pd(t(this.bee), this.tree, include.root = TRUE)
# })
# 
# PD <- do.call(rbind, PD)
# 
# spec.microbes <- cbind(spec.microbes, PD)
# 
# spec.net <- merge(spec.net, spec.microbes, all.x=TRUE)
# 
# ## check which individuals don't have microbe data
# drop.PD.NA <- unique(spec.net$UniqueID[spec.net$WeightsPar == 1 &
#                                          is.na(spec.net$PD)])
# 
# ## drop individuals that had parasite screen but not microbe
# ## filling in zeros for NAs in PD
# spec.net <- spec.net[!(spec.net$UniqueID %in% drop.PD.NA),] %>%
#   mutate(PD = ifelse(!is.na(PD), PD, 0))
# 
# 
# ##load tree from :
# ##Henriquez Piskulich, Patricia Andrea; Hugall, Andrew F.; Stuart-Fox, Devi (2023). A supermatrix phylogeny of the world’s bees (Hymenoptera: Anthophila) [Dataset]. Dryad. https://doi.org/10.5061/dryad.80gb5mkw1
# 
# phylo <- ape::read.tree("../../data/BEE_mat7_fulltree.nwk")
# 
# 
# ## i think we need to drop the tips that aren't in our dataset
# ## changing sp labels to match phylogeny
# spec.net$GenusSpecies <- gsub("Megachile gemula fulvogemula", "Megachile gemula", spec.net$GenusSpecies)
# spec.net$GenusSpecies <- gsub("Megachile melanophaea rohweri", "Megachile melanophaea", spec.net$GenusSpecies)
# 
# 
# ##clean up unwanted portion of labels
# pattern <- "(_n\\d+m\\d+_[A-Za-z0-9]+)?$"
# phylo$tip.label <- gsub(pattern, "", phylo$tip.label)
# 
# ## replace underscore with space
# phylo$tip.label <- gsub("_", " ", phylo$tip.label)
# 
# ## i think we need to drop the tips that aren't in our dataset
# ## changing sp labels to match phylogeny
# species_to_keep <- na.omit(unique(spec.net$GenusSpecies))
# 
# ## megachile comate and megachile subexilis are not in phylogeny so will drop these
# species_to_keep <- species_to_keep[!species_to_keep %in% c("Megachile comata", "Megachile subexilis")]
# 
# phylo_tips <- phylo$tip.label
# #only keep tips that match our species
# phylo <- ape::keep.tip(phylo, species_to_keep[species_to_keep %in% phylo_tips])
# 
# phylo_matrix <- ape::vcv.phylo(phylo)
# 
# 
# ## dropping species not in the phylogeny from the dataset
# ## megachile comate and megachile subexilis are not in phylogeny so will drop these
# spec.net <- spec.net[!spec.net$GenusSpecies %in% c("Megachile comata", "Megachile subexilis"),]
# # spec.net$GenusSpecies[spec.net$GenusSpecies == 'NA'] <- NA
# spec.net <- spec.net %>%
#   filter(is.na(GenusSpecies) == FALSE)
# 
# ## QUESTION: should there be NAs? not sure what check ids means
# ## check ids
# unique(spec.net$GenusSpecies[spec.net$Apidae == 1 &
#                                is.na(spec.net$MeanITD)])

#dropping genus species with NA -- 
# spec.net <- spec.net[!is.na(spec.net$GenusSpecies),]

## adding one to bombus abundance to make it non negative
# spec.net$Net_BombusAbundance <- spec.net$Net_BombusAbundance + 1
# spec.net$Net_HBAbundance <- spec.net$Net_HBAbundance + 1
# spec.net$Net_NonBombusHBAbundance <- spec.net$Net_NonBombusHBAbundance + 1


#squared transformed net bee abundance
#spec.net$Net_BeeAbundance <- (spec.net$Net_BeeAbundance)^2

## making site a factor so it works with gamma dist
#spec.net$Site <- as.factor(spec.net$Site)

##  center all of the x variables, need to use unique values to avoid
##  repetition by the number of specimens

#spec.net <- standardizeVars(spec.net, vars_yearsr, "YearSR")



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



## **********************************************************
## Bee abundance
## **********************************************************



## bee abund total
tot.bee.abund.vars <- c("MeanFloralAbundance",
                        "Year",
                        "SRDoy",
                        "I(SRDoy^2)",
                        "(1|Site)")

tot.bee.abund.x <- paste(tot.bee.abund.vars, collapse="+")
tot.bee.abund.y <- "BeeAbundance | weights(Weights)"
formula.tot.bee.abund <- as.formula(paste(tot.bee.abund.y, "~",tot.bee.abund.x))




#net bee abund
## bee abund total
net.bee.abund.vars <- c("MeanFloralAbundance",
                        "Year",
                        "SRDoy",
                        "I(SRDoy^2)",
                        "(1|Site)")

net.bee.abund.x <- paste(net.bee.abund.vars, collapse="+")
net.bee.abund.y <- "Net_BeeAbundance | weights(Weights)"
formula.net.bee.abund <- as.formula(paste(net.bee.abund.y, "~",net.bee.abund.x))





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

## bee div total
tot.bee.div.vars <- c("MeanFloralDiversity",
                      "Year",
                      "SRDoy",
                      "I(SRDoy^2)",
                      "Lat",
                      "(1|Site)")

tot.bee.div.x <- paste(tot.bee.div.vars, collapse="+")
tot.bee.div.y <- "BeeDiversity | weights(Weights)"
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
# bform.community <- bf.fabund +
#   bf.fdiv +
#   bf.tot.babund +
#   bf.tot.bdiv  +
#   set_rescor(FALSE)
# 
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
#                    "BeeDiversity", "Lat", #check this doesn't make VIF high
#                    "MeanFloralDiversity", "MeanITD", "Sociality", # if not at the genus level
#                    "(1|Site)", "rare.degree", "(1|gr(GenusSpecies, cov = phylo_matrix))") # add cov matrix for each genus
# # rare.degree
# # split genus into separate PD for apis bombus megachile
# 
# 
# microbe.x <- paste(microbe.vars, collapse="+")
# microbe.y <- "PD | weights(LogWeightsAbund)"
# formula.microbe <- as.formula(paste(microbe.y, "~",
#                                     microbe.x))
# 
# 
# bf.microbe <- bf(formula.microbe)

# #combine forms
# bform <- bf.fabund +
#   bf.fdiv +
#   bf.tot.babund +
#   bf.tot.bdiv  +
#   bf.microbe +
#   set_rescor(FALSE)
# 
# ## run model
# fit.microbe <- brm(bform , spec.net,
#                    cores=ncores,
#                    iter = 5000,
#                    chains =1,
#                    thin=1,
#                    init=0,
#                    open_progress = FALSE,
#                    control = list(adapt_delta = 0.99),
#                    save_pars = save_pars(all = TRUE),
#                    data2 = list(phylo_matrix=phylo_matrix))
# 
# write.ms.table(fit.microbe, "full_microbe")
# r2loo <- loo_R2(fit.microbe)
# r2 <- rstantools::bayes_R2(fit.microbe)
# save(fit.microbe, spec.net, r2, r2loo,
#      file="saved/fullMicrobeFit.Rdata")

## run bombus model
microbe.bombus.vars <- c("BeeAbundance",
                                   "BeeDiversity", "Lat", #check this doesn't make VIF high
                                   "MeanFloralDiversity", "MeanITD",  "rare.degree",
                                   "(1|Site)", "(1|gr(GenusSpecies, cov = phylo_matrix))") # add cov matrix for each genus
## NA check
 check_for_NA(microbe.bombus.vars)


microbe.bombus.x <- paste(microbe.bombus.vars, collapse="+")
microbe.bombus.y <- "PD | weights(LogWeightsAbund)"
formula.microbe.bombus <- as.formula(paste(microbe.bombus.y, "~",
                                     microbe.bombus.x))


bf.microbe.bombus <- bf(formula.microbe.bombus)

#combine forms
bform.bombus <- bf.fabund +
   bf.fdiv +
   bf.tot.babund +
   bf.tot.bdiv  +
   bf.microbe.bombus +
   set_rescor(FALSE)


fit.microbe.bombus <- brm(bform.bombus , spec.bombus,
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

## run apis model
# microbe.apis.vars <- c("BeeAbundance",
#                        "BeeDiversity", "Lat", #check this doesn't make VIF high
#                        "MeanFloralDiversity",# "MeanITD",
#                        "(1|Site)", "rare.degree"#, "(1|gr(GenusSpecies, cov = phylo_matrix))") # add cov matrix for each genus
# )
# 
# 
# microbe.apis.x <- paste(microbe.apis.vars, collapse="+")
# microbe.apis.y <- "PD | weights(LogWeightsAbund)"
# formula.microbe.apis <- as.formula(paste(microbe.apis.y, "~",
#                                          microbe.apis.x))
# 
# 
# bf.microbe.apis <- bf(formula.microbe.apis)
# 
# #combine forms
# bform.apis <- bf.fabund +
#   bf.fdiv +
#   bf.tot.babund +
#   bf.tot.bdiv  +
#   bf.microbe.apis +
#   set_rescor(FALSE)
# 
# fit.microbe.apis <- brm(bform.apis , spec.apis,
#                         cores=ncores,
#                         iter = 10000,
#                         chains =1,
#                         thin=1,
#                         init=0,
#                         open_progress = FALSE,
#                         control = list(adapt_delta = 0.99),
#                         save_pars = save_pars(all = TRUE))
# 
# write.ms.table(fit.microbe.apis, "apis_microbe")
# r2loo.apis <- loo_R2(fit.microbe.apis)
# r2.apis <- rstantools::bayes_R2(fit.microbe.apis)
# save(fit.microbe.apis, spec.apis, r2.apis, r2loo.apis,
#      file="saved/fullMicrobeApisFit.Rdata")

# ## run melissodes model
# microbe.melissodes.vars <- c("BeeAbundance",
#                          "BeeDiversity", "Lat", #check this doesn't make VIF high
#                          "MeanFloralDiversity", #"MeanITD",
#                          "(1|Site)", "rare.degree") #, "(1|gr(GenusSpecies, cov = phylo_matrix))") # add cov matrix for each genus
# 
# 
# 
# microbe.melissodes.x <- paste(microbe.melissodes.vars, collapse="+")
# microbe.melissodes.y <- "PD | weights(LogWeightsAbund)"
# formula.microbe.melissodes <- as.formula(paste(microbe.melissodes.y, "~",
#                                            microbe.melissodes.x))
# 
# 
# bf.microbe.melissodes <- bf(formula.microbe.melissodes)
# 
# #combine forms
# bform.melissodes <- bf.fabund +
#   bf.fdiv +
#   bf.tot.babund +
#   bf.tot.bdiv  +
#   bf.microbe.melissodes +
#   set_rescor(FALSE)
# 
# fit.microbe.melissodes <- brm(bform.melissodes , spec.melissodes,
#                           cores=ncores,
#                           iter = 10000,
#                           chains =1,
#                           thin=1,
#                           init=0,
#                           open_progress = FALSE,
#                           control = list(adapt_delta = 0.99),
#                           save_pars = save_pars(all = TRUE))
# 
# write.ms.table(fit.microbe.melissodes, "melissodes_microbe")
# r2loo.melissodes <- loo_R2(fit.microbe.melissodes)
# r2.melissodes <- rstantools::bayes_R2(fit.microbe.melissodes)
# save(fit.microbe.melissodes, spec.melissodes, r2.melissodes, r2loo.melissodes,
#      file="saved/fullMicrobeMelissodesFit.Rdata")



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