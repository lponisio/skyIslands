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
labs.itd <- (pretty(spec.sp$MeanITD, n=10))
axis.itd <-  standardize.axis(labs.itd,
                              spec.sp$MeanITD)


## load model output data
load(file="saved/fullMicrobeBombusFit.Rdata")
#load(file="saved/fullMicrobeApisFit.Rdata")
load(file="saved/fullMicrobeMelissodesFit.Rdata")


# We want the standarized data for the predictions (spec.data)
spec.bombus <- spec.net[spec.net$Genus == "Bombus",]
spec.melissodes <- spec.net[spec.net$Genus == "Melissodes",]
bombus.obligate <- spec.bombus[spec.bombus$WeightsObligateMicrobe == 1,]
data.site <- spec.net[spec.net$Weights == 1,]
bombus.transient <- spec.bombus[spec.bombus$WeightsTransientMicrobe == 1,]
melissodes.obligate <- spec.melissodes[spec.melissodes$WeightsObligateMicrobe == 1,]


## https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/


##CIs are 50-80-95

## **********************************************************
## Bombus model
## **********************************************************

## PD obligate ~ Bee Diversity

# Extract the data from conditional_effects
cond_effects_data <- conditional_effects(fit.microbe.bombus, effects = "BeeDiversity", resp = "PDobligate", plot = FALSE)
plot_data <- cond_effects_data$PDobligate.PDobligate_BeeDiversity


# Plot using ggplot2 for credible intervals with geom_ribbon
ob_beediv <- ggplot(plot_data, aes(x = BeeDiversity, y = estimate__)) +
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
  geom_point(data = bombus.obligate, aes(x = BeeDiversity, y = PD.obligate), 
             color = "black", alpha = 0.6) +
  # Labels and theme
  labs(x = "Bee Diversity (Untransformed)", y = "PD Obligate Estimate") +
  scale_x_continuous(breaks = axis.bee.div, labels = labs.bee.div) +
  theme_classic()

ob_beediv

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
  labs(x = "Diet breadth (rarefied degree; untransformed)", y = "PD Obligate Estimate") +
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
  labs(x = "Body Size (untransformed)", y = "PD Obligate Estimate") +
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
  labs(x = "Bee Abundance (Untransformed)", y = "PD Obligate Estimate") +
  scale_x_continuous(breaks = axis.bee.abund, labels = labs.bee.abund) +
  theme_classic()

ob_beeabund