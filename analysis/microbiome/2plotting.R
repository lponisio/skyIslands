## setwd('C:/Users/na_ma/Dropbox (University of Oregon)/Rotation/skyIslands')
## setwd('~/Dropbox (University of Oregon)/skyIslands')
setwd('~/Dropbox (University of Oregon)/skyIslands')
#setwd("/Volumes/bombus/rhayes/Dropbox (University of Oregon)/skyIslands")

setwd("analysis/microbiome/")

rm(list=ls())



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
make.plots = TRUE
source("src/init_microbe.R")
source("src/writeResultsTable.R")
source("src/runParasiteModels.R")
source("src/runPlotFreqModelDiagnostics.R")


## set to the number of cores you would like the models to run on
ncores <- 1

spec.net.orig <- spec.net

spec.bombus.orig <- spec.bombus
spec.apis.orig <- spec.apis
spec.melissodes.orig <- spec.melissodes

## load model output data
load(file="saved/fullMicrobeBombusFit.Rdata")
load(file="saved/fullMicrobeApisFit.Rdata")
load(file="saved/fullMicrobeMelissodesFit.Rdata")


##CIs are 50-80-95

## **********************************************************
## Bombus model
## **********************************************************

data.par <- spec.net[spec.net$BombusWeights != 0, ]
#
## PD ~ bee abundance
labs.bee.abund <- (pretty(c(0, spec.bombus.orig$BeeAbundance), n=8))
axis.bee.abund <-  standardize.axis(labs.bee.abund,
                                    spec.bombus.orig$BeeAbundance)

newdata.beeabund <- tidyr::crossing(BeeAbundance =
                                      seq(min(data.par$BeeAbundance),
                                          max(data.par$BeeAbundance),
                                          length.out=10),
                                    MeanFloralDiversity=mean(data.par$MeanFloralDiversity),
                                    Site=data.par$Site,
                                    rare.degree=mean(data.par$rare.degree),
                                    MeanITD=mean(data.par$MeanITD),
                                    Year = data.par$Year,
                                    DoyStart = data.par$DoyStart,
                                    Lat = mean(data.par$Lat),
                                    BeeDiversity = mean(data.par$BeeDiversity),
)

pred_beeabund <- fit.microbe.bombus %>%
  epred_draws(newdata = newdata.beeabund ,
              resp = "PD",
              allow_new_levels = TRUE)

## to see range of predicted values
# pred_beeabund %>%
#   group_by(BeeAbundance) %>%
#   summarise(mean(.epred))


bombus.abund.PD <- ggplot(pred_beeabund, aes(x = BeeAbundance, y =
                                                     .epred)) +
  stat_lineribbon(show.legend = FALSE) +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Bee Abundance", y = "Microbe \nPhylogenetic \nDistance",
       fill = "Credible Interval") +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16),
        legend.position = "none") +
  theme_classic() +
  geom_point(data=data.par,
             aes(y=PD, x=BeeAbundance), cex=2, alpha=0.5) +
  scale_x_continuous(
    breaks = axis.bee.abund,
    labels =  labs.bee.abund) + labs(tag='C.', y='Microbe \nPhylogenetic \nDistance')


bombus.abund.PD
panelC <- bombus.abund.PD

dir.create("figures")

ggsave(bombus.abund.PD, file="figures/bombusPD_beeAbund.pdf",
       height=4, width=5)

# # #################################
# ## PD ~ bee diversity
# labs.bee.div <- (pretty(c(0, spec.bombus.orig$BeeDiversity), n=8))
# axis.bee.div <-  standardize.axis(labs.bee.div,
#                                     spec.bombus.orig$BeeDiversity)
# 
# newdata.beediv <- tidyr::crossing(BeeDiversity =
#                                       seq(min(data.par$BeeDiversity),
#                                           max(data.par$BeeDiversity),
#                                           length.out=10),
#                                     MeanFloralDiversity=mean(data.par$MeanFloralDiversity),
#                                     Site=data.par$Site,
#                                     rare.degree=mean(data.par$rare.degree),
#                                     MeanITD=mean(data.par$MeanITD),
#                                     Year = data.par$Year,
#                                     DoyStart = data.par$DoyStart,
#                                     Lat = mean(data.par$Lat),
#                                     BeeAbundance = mean(data.par$BeeAbundance),
# )
# 
# pred_beediv <- fit.microbe.bombus %>%
#   epred_draws(newdata = newdata.beediv ,
#               resp = "PD",
#               allow_new_levels = TRUE)
# 
# 
# bombus.div.PD <- ggplot(pred_beediv, aes(x = BeeDiversity, y =
#                                                .epred)) +
#   stat_lineribbon(show.legend=FALSE) +
#   scale_fill_brewer(palette = "Oranges") +
#   labs(x = "Bee Diversity", y = "Microbe  \nPhylogenetic \nDistance",
#        fill = "Credible Interval") +
#   theme(legend.position = "none")  +
#   theme(axis.title.x = element_text(size=16),
#         axis.title.y = element_text(size=16),
#         text = element_text(size=16)) +
#   theme_classic() +
#   geom_point(data=data.par,
#              aes(y=PD, x=BeeDiversity), cex=2, alpha=0.5) +
#   scale_x_continuous(
#     breaks = axis.bee.div,
#     labels =  labs.bee.div) + labs(tag='C.', y='Microbe \nPhylogenetic \nDistance')
# 
# 
# bombus.div.PD
# #panelC <- bombus.div.PD
# 
# dir.create("figures")
# 
# ggsave(bombus.div.PD, file="figures/bombusPD_beeDiv.pdf",
#        height=4, width=5)
# 
# 
# # #################################
# ## PD ~ mean floral diversity
# labs.floral.div <- (pretty(c(0, spec.bombus.orig$MeanFloralDiversity), n=8))
# axis.floral.div <-  standardize.axis(labs.floral.div,
#                                   spec.bombus.orig$MeanFloralDiversity)
# 
# newdata.beediv <- tidyr::crossing(MeanFloralDiversity =
#                                     seq(min(data.par$MeanFloralDiversity),
#                                         max(data.par$MeanFloralDiversity),
#                                         length.out=10),
#                                   BeeDiversity=mean(data.par$BeeDiversity),
#                                   Site=data.par$Site,
#                                   rare.degree=mean(data.par$rare.degree),
#                                   MeanITD=mean(data.par$MeanITD),
#                                   Year = data.par$Year,
#                                   DoyStart = data.par$DoyStart,
#                                   Lat = mean(data.par$Lat),
#                                   BeeAbundance = mean(data.par$BeeAbundance),
# )
# 
# pred_beediv <- fit.microbe.bombus %>%
#   epred_draws(newdata = newdata.beediv ,
#               resp = "PD",
#               allow_new_levels = TRUE)
# 
# 
# floral.div.PD <- ggplot(pred_beediv, aes(x = MeanFloralDiversity, y =
#                                            .epred)) +
#   stat_lineribbon(show.legend=FALSE) +
#   scale_fill_brewer(palette = "Blues") +
#   labs(x = "Mean Floral Diversity", y = "Microbe Phylogenetic Diversity",
#        fill = "Credible Interval") +
#   theme(legend.position = "none")  +
#   theme(axis.title.x = element_text(size=16),
#         axis.title.y = element_text(size=16),
#         text = element_text(size=16)) +
#   theme_classic() +
#   geom_point(data=data.par,
#              aes(y=PD, x=MeanFloralDiversity), cex=2, alpha=0.5) +
#   scale_x_continuous(
#     breaks = axis.floral.div,
#     labels =  labs.floral.div) + labs(tag="A.", y='Microbe \nPhylogenetic \nDistance')
# 
# 
# floral.div.PD
# panelA <- floral.div.PD
# 
# dir.create("figures")
# 
# ggsave(floral.div.PD, file="figures/bombusPD_floralDiv.pdf",
#        height=4, width=5)
# 
# # #################################
# # ## PD ~ MeanITD
# labs.meanITD <- (pretty(c(0, spec.bombus.orig$MeanITD), n=8))
# axis.meanITD <-  standardize.axis(labs.meanITD,
#                                      spec.bombus.orig$MeanITD)
# 
# newdata.beediv <- tidyr::crossing(MeanITD =
#                                     seq(min(data.par$MeanITD),
#                                         max(data.par$MeanITD),
#                                         length.out=10),
#                                   BeeDiversity=mean(data.par$BeeDiversity),
#                                   Site=data.par$Site,
#                                   rare.degree=mean(data.par$rare.degree),
#                                   MeanFloralDiversity=mean(data.par$MeanFloralDiversity),
#                                   Year = data.par$Year,
#                                   DoyStart = data.par$DoyStart,
#                                   Lat = mean(data.par$Lat),
#                                   BeeAbundance = mean(data.par$BeeAbundance),
# )
# 
# pred_beediv <- fit.microbe.bombus %>%
#   epred_draws(newdata = newdata.beediv ,
#               resp = "PD",
#               allow_new_levels = TRUE)
# 
# 
# meanITD.PD <- ggplot(pred_beediv, aes(x = MeanITD, y =
#                                            .epred)) +
#   stat_lineribbon(show.legend=FALSE) +
#   scale_fill_brewer(palette = "Blues") +
#   labs(x = "Body Size (Mean ITD)", y = "Microbe Phylogenetic Diversity",
#        fill = "Credible Interval") +
#   theme(legend.position = "none")  +
#   theme(axis.title.x = element_text(size=16),
#         axis.title.y = element_text(size=16),
#         text = element_text(size=16)) +
#   theme_classic() +
#   geom_point(data=data.par,
#              aes(y=PD, x=MeanITD), cex=2, alpha=0.5) +
#   scale_x_continuous(
#     breaks = axis.meanITD,
#     labels =  labs.meanITD) + labs(tag='D.', y='Microbe \nPhylogenetic \nDistance')
# 
# 
# meanITD.PD
# panelD <- meanITD.PD
# 
# dir.create("figures")
# 
# ggsave(meanITD.PD, file="figures/bombusPD_meanITD.pdf",
#        height=4, width=5)
# 
# # #################################
# # ## PD ~ rare.degree
# labs.rare.degree <- (pretty(c(0, spec.bombus.orig$rare.degree), n=8))
# axis.rare.degree <-  standardize.axis(labs.rare.degree,
#                                   spec.bombus.orig$rare.degree)
# 
# newdata.beediv <- tidyr::crossing(rare.degree =
#                                     seq(min(data.par$rare.degree),
#                                         max(data.par$rare.degree),
#                                         length.out=10),
#                                   BeeDiversity=mean(data.par$BeeDiversity),
#                                   Site=data.par$Site,
#                                   MeanITD=mean(data.par$MeanITD),
#                                   MeanFloralDiversity=mean(data.par$MeanFloralDiversity),
#                                   Year = data.par$Year,
#                                   DoyStart = data.par$DoyStart,
#                                   Lat = mean(data.par$Lat),
#                                   BeeAbundance = mean(data.par$BeeAbundance),
# )
# 
# pred_beediv <- fit.microbe.bombus %>%
#   epred_draws(newdata = newdata.beediv ,
#               resp = "PD",
#               allow_new_levels = TRUE)
# 
# 
# rare.degree.PD <- ggplot(pred_beediv, aes(x = rare.degree, y =
#                                         .epred)) +
#   stat_lineribbon(show.legend=FALSE) +
#   scale_fill_brewer(palette = "Blues") +
#   labs(x = "Rarefied Degree", y = "Microbe Phylogenetic Diversity",
#        fill = "Credible Interval") +
#   theme(legend.position = "none")  +
#   theme(axis.title.x = element_text(size=16),
#         axis.title.y = element_text(size=16),
#         text = element_text(size=16)) +
#   theme_classic() +
#   geom_point(data=data.par,
#              aes(y=PD, x=rare.degree), cex=2, alpha=0.5) +
#   scale_x_continuous(
#     breaks = axis.rare.degree,
#     labels =  labs.rare.degree) + labs(tag='B.', y='Microbe \nPhylogenetic \nDistance')
# 
# 
# rare.degree.PD
# panelB <- rare.degree.PD
# 
# dir.create("figures")
# 
# ggsave(rare.degree.PD, file="figures/bombusPD_rare.degree.pdf",
#        height=4, width=5)
# 
# #make paneled fig
# grid.arrange(panelA, 
#              panelB,
#              panelC, 
#              ncol=1)



## **********************************************************
## Apis model
## **********************************************************
# 
data.par <- spec.apis[spec.apis$ApisWeights != 0, ]

## PD ~ bee abundance
labs.bee.abund <- (pretty(c(0, spec.apis.orig$BeeAbundance), n=8))
axis.bee.abund <-  standardize.axis(labs.bee.abund,
                                    spec.apis.orig$BeeAbundance)

newdata.beeabund <- tidyr::crossing(BeeAbundance =
                                      seq(min(data.par$BeeAbundance),
                                          max(data.par$BeeAbundance),
                                          length.out=10),
                                    MeanFloralDiversity=mean(data.par$MeanFloralDiversity),
                                    Site=data.par$Site,
                                    rare.degree=mean(data.par$rare.degree),
                                    #MeanITD=mean(data.par$MeanITD),
                                    Year = data.par$Year,
                                    DoyStart = data.par$DoyStart,
                                    Lat = mean(data.par$Lat),
                                    BeeDiversity = mean(data.par$BeeDiversity),
)

pred_beeabund <- fit.microbe.apis %>%
  epred_draws(newdata = newdata.beeabund ,
              resp = "PD",
              allow_new_levels = TRUE)

## to see range of predicted values
#pred_beeabund %>%
#  group_by(BeeAbundance) %>%
#  summarise(mean(.epred))


apis.abund.PD <- ggplot(pred_beeabund, aes(x = BeeAbundance, y =
                                               .epred)) +
  stat_lineribbon(show.legend=FALSE) +
  scale_fill_brewer(palette = "Oranges") +
  labs(x = "Bee Abundance", y = "Microbe \nPhylogenetic \nDistance",
       fill = "Credible Interval", tag="C.") +
  theme(legend.position = "none")  +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_classic() +
  geom_point(data=data.par,
             aes(y=PD, x=BeeAbundance), cex=2, alpha=0.5) +
  scale_x_continuous(
    breaks = axis.bee.abund,
    labels =  labs.bee.abund)


apis.abund.PD
panelC <- apis.abund.PD

dir.create("figures")

ggsave(apis.abund.PD, file="figures/apisPD_beeAbund.pdf",
       height=4, width=5)

#################################
## PD ~ bee diversity
labs.bee.div <- (pretty(c(0, spec.apis.orig$BeeDiversity), n=8))
axis.bee.div <-  standardize.axis(labs.bee.div,
                                  spec.apis.orig$BeeDiversity)

newdata.beediv <- tidyr::crossing(BeeDiversity =
                                    seq(min(data.par$BeeDiversity),
                                        max(data.par$BeeDiversity),
                                        length.out=10),
                                  MeanFloralDiversity=mean(data.par$MeanFloralDiversity),
                                  Site=data.par$Site,
                                  rare.degree=mean(data.par$rare.degree),
                                  #MeanITD=mean(data.par$MeanITD),
                                  Year = data.par$Year,
                                  DoyStart = data.par$DoyStart,
                                  Lat = mean(data.par$Lat),
                                  BeeAbundance = mean(data.par$BeeAbundance),
)

pred_beediv <- fit.microbe.apis %>%
  epred_draws(newdata = newdata.beediv ,
              resp = "PD",
              allow_new_levels = TRUE)


apis.div.PD <- ggplot(pred_beediv, aes(x = BeeDiversity, y =
                                           .epred)) +
  stat_lineribbon(show.legend=FALSE) +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Bee Diversity", y = "Microbe \nPhylogenetic \nDistance",
       fill = "Credible Interval", tag="A.") +
  theme(legend.position = "none")  +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_classic() +
  geom_point(data=data.par,
             aes(y=PD, x=BeeDiversity), cex=2, alpha=0.5) +
  scale_x_continuous(
    breaks = axis.bee.div,
    labels =  labs.bee.div)


apis.div.PD
panelA <- apis.div.PD

dir.create("figures")

ggsave(apis.div.PD, file="figures/apisPD_beeDiv.pdf",
       height=4, width=5)


#################################
## PD ~ mean floral diversity
labs.floral.div <- (pretty(c(0, spec.apis.orig$MeanFloralDiversity), n=8))
axis.floral.div <-  standardize.axis(labs.floral.div,
                                     spec.apis.orig$MeanFloralDiversity)

newdata.flrdiv <- tidyr::crossing(MeanFloralDiversity =
                                    seq(min(data.par$MeanFloralDiversity),
                                        max(data.par$MeanFloralDiversity),
                                        length.out=10),
                                  BeeDiversity=mean(data.par$BeeDiversity),
                                  Site=data.par$Site,
                                  rare.degree=mean(data.par$rare.degree),
                                  #MeanITD=mean(data.par$MeanITD),
                                  Year = data.par$Year,
                                  DoyStart = data.par$DoyStart,
                                  Lat = mean(data.par$Lat),
                                  BeeAbundance = mean(data.par$BeeAbundance),
)

pred_flrdiv <- fit.microbe.apis %>%
  epred_draws(newdata = newdata.flrdiv ,
              resp = "PD",
              allow_new_levels = TRUE)


floral.div.PD <- ggplot(pred_flrdiv, aes(x = MeanFloralDiversity, y =
                                           .epred)) +
  stat_lineribbon(show.legend=FALSE) +
  scale_fill_brewer(palette = "Oranges") +
  labs(x = "Mean Floral Diversity", y = "Microbe \nPhylogenetic \nDistance",
       fill = "Credible Interval", tag="D.") +
  theme(legend.position = "none")  +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_classic() +
  geom_point(data=data.par,
             aes(y=PD, x=MeanFloralDiversity), cex=2, alpha=0.5) +
  scale_x_continuous(
    breaks = axis.floral.div,
    labels =  labs.floral.div)


floral.div.PD
panelD <- floral.div.PD

dir.create("figures")

ggsave(floral.div.PD, file="figures/apisPD_floralDiv.pdf",
       height=4, width=5)


#################################
## PD ~ rare.degree
labs.rare.degree <- (pretty(c(0, spec.apis.orig$rare.degree), n=8))
axis.rare.degree <-  standardize.axis(labs.rare.degree,
                                      spec.apis.orig$rare.degree)

newdata.beediv <- tidyr::crossing(rare.degree =
                                    seq(min(data.par$rare.degree),
                                        max(data.par$rare.degree),
                                        length.out=10),
                                  BeeDiversity=mean(data.par$BeeDiversity),
                                  Site=data.par$Site,
                                  #MeanITD=mean(data.par$MeanITD),
                                  MeanFloralDiversity=mean(data.par$MeanFloralDiversity),
                                  Year = data.par$Year,
                                  DoyStart = data.par$DoyStart,
                                  Lat = mean(data.par$Lat),
                                  BeeAbundance = mean(data.par$BeeAbundance),
)

pred_beediv <- fit.microbe.apis %>%
  epred_draws(newdata = newdata.beediv ,
              resp = "PD",
              allow_new_levels = TRUE)


rare.degree.PD <- ggplot(pred_beediv, aes(x = rare.degree, y =
                                            .epred)) +
  stat_lineribbon(show.legend=FALSE) +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Rarified Degree", y = "Microbe \nPhylogenetic \nDistance",
       fill = "Credible Interval", tag="B.") +
  theme(legend.position = "none")  +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_classic() +
  geom_point(data=data.par,
             aes(y=PD, x=rare.degree), cex=2, alpha=0.5) +
  scale_x_continuous(
    breaks = axis.rare.degree,
    labels =  labs.rare.degree)


rare.degree.PD
panelB <- rare.degree.PD

dir.create("figures")

ggsave(rare.degree.PD, file="figures/apisPD_rare.degree.pdf",
       height=4, width=5)

#################################
## PD ~ latitude
labs.lat <- (pretty(c(0, spec.apis.orig$Lat), n=8))
axis.lat <-  standardize.axis(labs.lat,
                                      spec.apis.orig$Lat)

newdata.beediv <- tidyr::crossing(Lat =
                                    seq(min(data.par$Lat),
                                        max(data.par$Lat),
                                        length.out=10),
                                  BeeDiversity=mean(data.par$BeeDiversity),
                                  Site=data.par$Site,
                                  rare.degree=mean(data.par$rare.degree),
                                  MeanFloralDiversity=mean(data.par$MeanFloralDiversity),
                                  Year = data.par$Year,
                                  DoyStart = data.par$DoyStart,
                                  BeeAbundance = mean(data.par$BeeAbundance),
)

pred_beediv <- fit.microbe.apis %>%
  epred_draws(newdata = newdata.beediv ,
              resp = "PD",
              allow_new_levels = TRUE)


lat.PD <- ggplot(pred_beediv, aes(x = Lat, y =
                                            .epred)) +
  stat_lineribbon(show.legend=FALSE) +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Latitude", y = "Microbe \nPhylogenetic \nDistance",
       fill = "Credible Interval", tag="E.") +
  theme(legend.position = "none")  +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_classic() +
  geom_point(data=data.par,
             aes(y=PD, x=Lat), cex=2, alpha=0.5) +
  scale_x_continuous(
    breaks = axis.lat,
    labels =  labs.lat)


lat.PD
panelE <- rare.lat.PD

dir.create("figures")

ggsave(rare.degree.PD, file="figures/apisPD_Lat.pdf",
       height=4, width=5)


##make paneled fig
grid.arrange(panelA, panelB,
             panelC, panelD, panelE, ncol=2)


## Melissodes model

# data.par <- spec.net[spec.net$MelissodesWeights != 0, ]
# 
# ## PD ~ bee abundance
# labs.bee.abund <- (pretty(c(0, spec.melissodes.orig$BeeAbundance), n=8))
# axis.bee.abund <-  standardize.axis(labs.bee.abund,
#                                     spec.melissodes.orig$BeeAbundance)
# 
# newdata.beeabund <- tidyr::crossing(BeeAbundance =
#                                       seq(min(data.par$BeeAbundance),
#                                           max(data.par$BeeAbundance),
#                                           length.out=10),
#                                     MeanFloralDiversity=mean(data.par$MeanFloralDiversity),
#                                     Site=data.par$Site,
#                                     rare.degree=mean(data.par$rare.degree),
#                                     #MeanITD=mean(data.par$MeanITD),
#                                     Year = data.par$Year,
#                                     DoyStart = data.par$DoyStart,
#                                     Lat = mean(data.par$Lat),
#                                     BeeDiversity = mean(data.par$BeeDiversity),
# )
# 
# pred_beeabund <- fit.microbe.melissodes %>%
#   epred_draws(newdata = newdata.beeabund ,
#               resp = "PD",
#               allow_new_levels = TRUE)
# 
# ## to see range of predicted values
# #pred_beeabund %>%
# #  group_by(BeeAbundance) %>%
# #  summarise(mean(.epred))
# 
# 
# melissodes.abund.PD <- ggplot(pred_beeabund, aes(x = BeeAbundance, y =
#                                              .epred)) +
#   stat_lineribbon(show.legend=FALSE) +
#   scale_fill_brewer(palette = "Oranges") +
#   labs(x = "Bee Abundance", y = "Microbe \nPhylogenetic \nDistance",
#        fill = "Credible Interval") +
#   theme(legend.position = "none")  +
#   theme(axis.title.x = element_text(size=16),
#         axis.title.y = element_text(size=16),
#         text = element_text(size=16)) +
#   theme_classic() +
#   geom_point(data=data.par,
#              aes(y=PD, x=BeeAbundance), cex=2, alpha=0.5) +
#   scale_x_continuous(
#     breaks = axis.bee.abund,
#     labels =  labs.bee.abund)
# 
# 
# melissodes.abund.PD
# 
# dir.create("figures")
# 
# ggsave(melissodes.abund.PD, file="figures/melissodesPD_beeAbund.pdf",
#        height=4, width=5)
# 
