rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path
setwd(local.path)
setwd("skyIslands/analysis/microbiome/")


library(lme4)
#library(glmmADMB)
library(R2admb)
library(shinystan)
library(picante)
library(plyr)
library(bayesplot)
library(pscl)
library(brms)
library(performance)
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
vars_sp <- "MeanITD.x"
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

## ***********************************************************************
## plotting, unscaling labs
## ***********************************************************************
## lat

spec.uni <- spec.orig[spec.orig$Weights ==1,]


## Bombus mods
##  PDobligate ~ BeeDiversity
##  PDobligate ~ rare.degree
##  PDtransientlog ~ MeanITD

## TODO fix labels for these
## TODO add code from untitled 1
## TODO add melissodes plots
## TODO add apis plots when done fixing mods
## TODO save out plots
## TODO prep for tomorrow meeting

## TODO ask difference between spec.uni and spec.sp
## bee diversity
labs.bee.abund <- (pretty(c(spec.uni$BeeAbundance), n=8))
axis.bee.abund <-  standardize.axis(labs.bee.abund,
                                  spec.uni$BeeAbundance)

## bee diversity
labs.bee.div <- (pretty(c(spec.uni$BeeDiversity), n=8))
axis.bee.div <-  standardize.axis(labs.bee.div,
                                  spec.uni$BeeDiversity)

spec.sp <-  spec.orig[spec.orig$WeightsObligateMicrobe == 1 & spec.orig$Genus == "Bombus",]

## rare.degree
labs.degree <- (pretty(spec.sp$rare.degree, n=10))
axis.degree <-  standardize.axis(labs.degree,
                                 spec.sp$rare.degree)

spec.sp <-  spec.orig[spec.orig$WeightsTransientMicrobe == 1 & spec.orig$Genus == "Bombus",]
## meanITD
labs.itd <- (pretty(spec.sp$MeanITD.x, n=4))
axis.itd <-  standardize.axis(labs.itd,
                              spec.sp$MeanITD.x)


## load model output data
load(file="saved/fullMicrobeBombusFit.Rdata")

#when model finishes change over
#load(file="saved/fullMicrobeBombusFit_2.Rdata")
load(file="saved/fullMicrobeMelissodesFit_2.Rdata")

#load(file="saved/fullMicrobeApisFit.Rdata")


# We want the standarized data for the predictions (spec.data)
spec.bombus <- spec.net[spec.net$Genus == "Bombus",]
bombus.obligate <- spec.bombus[spec.bombus$WeightsObligateMicrobe == 1,]
bombus.transient <- spec.bombus[spec.bombus$WeightsTransientMicrobe == 1,]

spec.melissodes <- spec.net[spec.net$Genus == "Melissodes",]
melissodes.obligate <- spec.melissodes[spec.melissodes$WeightsObligateMicrobe == 1,]
melissodes.transient <- spec.melissodes[spec.melissodes$WeightsTransientMicrobe == 1,]


## https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/


##CIs are 50-80-95

## **********************************************************
## Bombus model
## **********************************************************


## TODO: check model results for bombus with log transformed PD to simplify plots, right now just going ahead with unmatched y axis


plot_model_condeff_compare <- function(model.a,
                                       model.b,
                                       this.effect,
                                       this.resp.a,
                                       this.resp.b,
                                       ){

## PD obligate ~ Bee Diversity

# Extract the data from conditional_effects
cond_effects_data_a <- conditional_effects(fit.microbe.bombus, effects = "BeeDiversity", resp = "PDobligate", plot = FALSE)
plot_data_a <- cond_effects_data_a$PDobligate.PDobligate_BeeDiversity

# Extract the data from conditional_effects
cond_effects_data_b <- conditional_effects(fit.microbe.melissodes, effects = "BeeDiversity", resp = "PDobligatelog", plot = FALSE)
plot_data_b <- cond_effects_data_b$PDobligatelog.PDobligatelog_BeeDiversity

# Plot using ggplot2 for credible intervals with geom_ribbon
ob_beediv <- ggplot(plot_data_a, aes(x = BeeDiversity, y = estimate__)) +
  # Add ribbons for the 95%, 80%, and 50% credible intervals
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill='navy', color = "navy", linetype='solid') +
  geom_ribbon(aes(ymin = lower__ + 0.1 * (upper__ - lower__), 
                  ymax = upper__ - 0.1 * (upper__ - lower__)), 
              alpha = 0.3, fill='navy', color = "navy", linetype='dashed') +
  geom_ribbon(aes(ymin = lower__ + 0.25 * (upper__ - lower__), 
                  ymax = upper__ - 0.25 * (upper__ - lower__)), 
              alpha = 0.4, fill='navy', color = "navy", linetype='dotted') +
  # Add ribbons for the 95%, 80%, and 50% credible intervals
  geom_ribbon(data=plot_data_b, aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill='orange', color = "orange", linetype='solid') +
  geom_ribbon(data=plot_data_b, aes(ymin = lower__ + 0.1 * (upper__ - lower__), 
                  ymax = upper__ - 0.1 * (upper__ - lower__)), 
              alpha = 0.3, fill='orange', color = "orange", linetype='dashed') +
  geom_ribbon(data=plot_data_b, aes(ymin = lower__ + 0.25 * (upper__ - lower__), 
                  ymax = upper__ - 0.25 * (upper__ - lower__)), 
              alpha = 0.4, fill='orange', color = "orange", linetype='dotted') +
  # Add line for the estimates
  geom_line(color = "navy") +
  # Add points for original data
  geom_point(data = bombus.obligate, aes(x = BeeDiversity, y = PD.obligate), 
             color = "navy", alpha = 0.6) +
  # Add line for the estimates
  geom_line(data=plot_data_b, aes(x = BeeDiversity, y = PD.obligate.log), color = "orange") +
  # Add points for original data
  geom_point(data = melissodes.obligate, aes(x = BeeDiversity, y = PD.obligate.log), 
             color = "orange", alpha = 0.6) +
  # Labels and theme
  labs(x = "Bee Diversity (Untransformed)", y = "Obligate Microbe PD") +
  scale_x_continuous(breaks = axis.bee.div, labels = labs.bee.div) +
  theme_classic()

ob_beediv
}
############################
## PDobligate ~ rare.degree

# Extract the data from conditional_effects
cond_effects_data <- conditional_effects(fit.microbe.bombus, effects = "rare.degree", resp = "PDobligate", plot = FALSE)
plot_data <- cond_effects_data$PDobligate.PDobligate_rare.degree


# Plot using ggplot2 for credible intervals with geom_ribbon
ob_degree <- ggplot(plot_data, aes(x = rare.degree, y = estimate__)) +
  # Add ribbons for the 95%, 80%, and 50% credible intervals
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill = "darkgreen") +
  geom_ribbon(aes(ymin = lower__ + 0.1 * (upper__ - lower__), 
                  ymax = upper__ - 0.1 * (upper__ - lower__)), 
              alpha = 0.3, fill = "darkgreen") +
  geom_ribbon(aes(ymin = lower__ + 0.25 * (upper__ - lower__), 
                  ymax = upper__ - 0.25 * (upper__ - lower__)), 
              alpha = 0.4, fill = "darkgreen") +
  # Add line for the estimates
  geom_line(color = "black") +
  # Add points for original data
  geom_point(data = bombus.obligate, aes(x = rare.degree, y = PD.obligate), 
             color = "black", alpha = 0.6) +
  # Labels and theme
  labs(x = "Diet breadth (rarefied degree; untransformed)", y = "Obligate Microbe PD") +
  scale_x_continuous(breaks = axis.degree, labels = labs.degree) +
  theme_classic()

ob_degree

############################
## PDtransientlog ~ MeanITD

# Extract the data from conditional_effects
cond_effects_data <- conditional_effects(fit.microbe.bombus, effects = "MeanITD", resp = "PDtransientlog", plot = FALSE)
plot_data <- cond_effects_data$PDtransientlog.PDtransientlog_MeanITD


# Plot using ggplot2 for credible intervals with geom_ribbon
trans_itd <- ggplot(plot_data, aes(x = MeanITD, y = estimate__)) +
  # Add ribbons for the 95%, 80%, and 50% credible intervals
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill = "orange") +
  geom_ribbon(aes(ymin = lower__ + 0.1 * (upper__ - lower__), 
                  ymax = upper__ - 0.1 * (upper__ - lower__)), 
              alpha = 0.3, fill = "orange") +
  geom_ribbon(aes(ymin = lower__ + 0.25 * (upper__ - lower__), 
                  ymax = upper__ - 0.25 * (upper__ - lower__)), 
              alpha = 0.4, fill = "orange") +
  # Add line for the estimates
  geom_line(color = "black") +
  # Add points for original data
  geom_point(data = bombus.transient, aes(x = MeanITD, y = PD.transient.log), 
             color = "black", alpha = 0.6) +
  # Labels and theme
  labs(x = "Body Size (untransformed)", y = "Facultative Microbe PD (log transformed)") +
  scale_x_continuous(breaks = axis.itd, labels = labs.itd) +
  theme_classic()

trans_itd

## **********************************************************
## Melissodes model
## **********************************************************

## PD obligate ~ Bee Diversity

# Extract the data from conditional_effects
cond_effects_data <- conditional_effects(fit.microbe.melissodes, effects = "BeeAbundance", resp = "PDobligatelog", plot = FALSE)
plot_data <- cond_effects_data$PDobligatelog.PDobligatelog_BeeAbundance


# Plot using ggplot2 for credible intervals with geom_ribbon
ob_beeabund <- ggplot(plot_data, aes(x = BeeAbundance, y = estimate__)) +
  # Add ribbons for the 95%, 80%, and 50% credible intervals
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill = "darkgreen") +
  geom_ribbon(aes(ymin = lower__ + 0.1 * (upper__ - lower__), 
                  ymax = upper__ - 0.1 * (upper__ - lower__)), 
              alpha = 0.3, fill = "darkgreen") +
  geom_ribbon(aes(ymin = lower__ + 0.25 * (upper__ - lower__), 
                  ymax = upper__ - 0.25 * (upper__ - lower__)), 
              alpha = 0.4, fill = "darkgreen") +
  # Add line for the estimates
  geom_line(color = "black") +
  # Add points for original data
  geom_point(data = melissodes.obligate, aes(x = BeeAbundance, y = PD.obligate.log), 
             color = "black", alpha = 0.6) +
  # Labels and theme
  labs(x = "Bee Abundance (Untransformed)", y = "Obligate Microbe PD (log transformed)") +
  scale_x_continuous(breaks = axis.bee.abund, labels = labs.bee.abund) +
  theme_classic()

ob_beeabund