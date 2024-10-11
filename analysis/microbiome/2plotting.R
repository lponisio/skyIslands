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
## bee abund
## TODO unlog transform the data for axes
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
load(file="saved/fullMicrobeBombusFit_2.Rdata")

#when model finishes change over
#load(file="saved/fullMicrobeBombusFit_2.Rdata")
load(file="saved/fullMicrobeMelissodesFit_2.Rdata")

#load(file="saved/fullMicrobeApisFit.Rdata")


# We want the standarized data for the predictions (spec.data)
spec.bombus <- spec.net[spec.net$Genus == "Bombus",]
bombus.obligate <- spec.bombus[spec.bombus$WeightsObligateMicrobe == 1,]
bombus.obligate$PDobligate <- bombus.obligate$PD.obligate
bombus.transient <- spec.bombus[spec.bombus$WeightsTransientMicrobe == 1,]
bombus.transient$PDtransient <- bombus.transient$PD.transient
bombus.transient$PDtransientlog <- bombus.transient$PD.transient.log

spec.melissodes <- spec.net[spec.net$Genus == "Melissodes",]
melissodes.obligate <- spec.melissodes[spec.melissodes$WeightsObligateMicrobe == 1,]
melissodes.obligate$PDobligate <- melissodes.obligate$PD.obligate
melissodes.obligate$PDobligatelog <- melissodes.obligate$PD.obligate.log
melissodes.transient <- spec.melissodes[spec.melissodes$WeightsTransientMicrobe == 1,]
melissodes.transient$PDtransient <- melissodes.transient$PD.transient
melissodes.transient$PDtransientlog <- melissodes.transient$PD.transient.log

## https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/


##CIs are 50-80-95

## **********************************************************
## Bombus model
## **********************************************************


## TODO: check model results for bombus with log transformed PD to simplify plots, right now just going ahead with unmatched y axis
source("src/plotting_functions.R")



obligate.bee.div.plot <- plot_model_condeff_compare(model.a=fit.microbe.bombus,
                                       model.b=fit.microbe.melissodes,
                                       this.effect='BeeDiversity',
                                       this.resp.a='PDobligate', ## TODO: potentially update to both be PD obligate log?
                                       this.resp.b='PDobligatelog',
                                       point.data.a=bombus.obligate,
                                       point.data.b=melissodes.obligate,
                                       axis.breaks=axis.bee.div,
                                       axis.labs=labs.bee.div,
                                       xlabel="Bee Diversity",
                                       ylabel="Obligate Microbe PD",
                                       mod1color='navy',
                                       mod2color='coral',
                                       fill.a=TRUE,
                                       fill.b=FALSE
                                       )
obligate.bee.div.plot 
############################
## PDobligate ~ rare.degree
  
  
obligate.rare.degree.plot <- plot_model_condeff_compare(model.a=fit.microbe.bombus,
                                         model.b=fit.microbe.melissodes,
                                         this.effect='rare.degree',
                                         this.resp.a='PDobligate', ## TODO: potentially update to both be PD obligate log?
                                         this.resp.b='PDobligatelog',
                                         point.data.a=bombus.obligate,
                                         point.data.b=melissodes.obligate,
                                         axis.breaks=axis.degree,
                                         axis.labs=labs.degree,
                                         xlabel="Diet Breadth",
                                         ylabel="Obligate Microbe PD",
                                         mod1color='navy',
                                         mod2color='coral',
                                         fill.a=TRUE,
                                         fill.b=FALSE
                                         )
obligate.rare.degree.plot 

############################
## PDtransientlog ~ MeanITD
## can only make for Bombus since melissodes doesn't have mean ITD

# obligate.ITD.plot <- plot_model_condeff_compare(model.a=fit.microbe.bombus,
#                                                         model.b=fit.microbe.melissodes,
#                                                         this.effect='MeanITD.x',
#                                                         this.resp.a="PDtransientlog", ## TODO: potentially update to both be PD obligate log?
#                                                         this.resp.b="PDtransientlog",
#                                                         point.data.a=bombus.transient,
#                                                         point.data.b=melissodes.transient,
#                                                         axis.breaks=axis.itd,
#                                                         axis.labs=labs.itd,
#                                                         xlabel="Body Size",
#                                                         ylabel="Transient Microbe PD (log)",
#                                                         mod1color='navy',
#                                                         mod2color='coral'
# )
# obligate.ITD.plot 
# Extract the data from conditional_effects
cond_effects_data <- conditional_effects(fit.microbe.bombus, effects = "MeanITD.x", resp = "PDtransientlog", plot = FALSE)
plot_data <- cond_effects_data$PDtransientlog.PDtransientlog_MeanITD.x


# Plot using ggplot2 for credible intervals with geom_ribbon
trans_itd <- ggplot(plot_data, aes(x = MeanITD.x, y = estimate__)) +
  # Add ribbons for the 95%, 80%, and 50% credible intervals
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill = "navy", linetype='dotted') +
  geom_ribbon(aes(ymin = lower__ + 0.1 * (upper__ - lower__), 
                  ymax = upper__ - 0.1 * (upper__ - lower__)), 
              alpha = 0.3, fill = "navy", linetype='dashed') +
  geom_ribbon(aes(ymin = lower__ + 0.25 * (upper__ - lower__), 
                  ymax = upper__ - 0.25 * (upper__ - lower__)), 
              alpha = 0.4, fill = "navy", linetype='solid') +
  # Add line for the estimates
  geom_line(color = "black") +
  # Add points for original data
  geom_point(data = bombus.transient, aes(x = MeanITD.x, y = PD.transient.log), 
             color = "navy", alpha = 0.6, size) +
  # Labels and theme
  labs(x = "Body Size", y = "Facultative Microbe PD") +
  scale_x_continuous(breaks = axis.itd, labels = labs.itd) +
  theme_classic()

trans_itd

## **********************************************************
## Melissodes model
## **********************************************************

## PD obligate log ~ Bee abundance
obligate.abund.plot <- plot_model_condeff_compare(model.a=fit.microbe.bombus,
                                                model.b=fit.microbe.melissodes,
                                                this.effect='BeeAbundance',
                                                this.resp.a="PDobligate", ## TODO: potentially update to both be PD obligate log?
                                                this.resp.b="PDobligatelog",
                                                point.data.a=bombus.obligate,
                                                point.data.b=melissodes.obligate,
                                                axis.breaks=axis.bee.abund,
                                                axis.labs=labs.bee.abund,
                                                xlabel="Bee Abundance",
                                                ylabel="Obligate Microbe PD",
                                                mod1color='navy',
                                                mod2color='coral',
                                                fill.a=FALSE,
                                                fill.b=TRUE
)
obligate.abund.plot


## TODO: make grid plots and save out to final


