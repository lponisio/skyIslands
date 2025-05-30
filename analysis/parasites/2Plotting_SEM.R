rm(list=ls())
## setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')
 setwd("C:/Users/na_ma/Dropbox (University of Oregon)/skyIslands")
##setwd('~/Dropbox (University of Oregon)/skyislands')
## Script for plotting all of the important explanatory variables.
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(brms)
library(tidybayes)
setwd("analysis/parasites")
load(file="saved/spec_weights.Rdata")
source("src/misc.R")
source("src/ggplotThemes.R")

## ***********************************************************************
## scaling/unscaling labs
## ***********************************************************************
# original data, subsetted to unique values
spec.uni.orig <- spec.orig[spec.orig$Weights ==1,]
## scaled data subsetted to unique values
spec.uni <- spec.net[spec.net$Weights ==1,]

## use unscaled data to have nice axis labels, convert to scaled for
## the axes

## lat (logged)
labs.lat.x <- pretty(c(spec.uni.orig$Lat),
                      n=10)
axis.lat.x <-  standardize.axis(labs.lat.x, spec.uni.orig$Lat)

## bloom abundance (not logged)
labs.bloom.abund <- (pretty(c(spec.uni.orig$MeanFloralAbundance), n=6))
axis.bloom.abund <-  standardize.axis(labs.bloom.abund,
                                      spec.uni.orig$MeanFloralAbundance)
## flower div (not logged)
labs.flower.div <- (pretty(spec.uni.orig$MeanFloralDiversity, n=5))
axis.flower.div <-  standardize.axis(labs.flower.div,
                                     spec.uni.orig$MeanFloralDiversity)
## HB abund (logged + 1)
labs.HB.abund <- (pretty(c(spec.uni.orig$Net_HBAbundance), n=5))
axis.HB.abund <-  standardize.axis(labs.HB.abund, spec.uni.orig$Net_HBAbundance)
## bombus abund (logged + 1)
labs.bombus.abund <- (pretty(c(spec.uni.orig$Net_BombusAbundance), n=5))
axis.bombus.abund <-  standardize.axis(labs.bombus.abund, spec.uni.orig$Net_BombusAbundance)
## all bee abund (logged)
labs.bee.abund <- (pretty(c(spec.uni.orig$Net_BeeAbundance), n=5))
axis.bee.abund <-  standardize.axis(labs.bee.abund, spec.uni.orig$Net_BeeAbundance)
## bee diversity (not logged)
labs.bee.div <- (pretty(c(spec.uni.orig$Net_BeeDiversity), n=5))
axis.bee.div <-  standardize.axis(labs.bee.div,
                                  spec.uni.orig$Net_BeeDiversity)

## use all the species data or just bombus? 
bombus.par <- spec.orig[spec.orig$WeightsPar==1 & spec.orig$Genus == "Bombus", ]
## rare.degree (logged)
labs.degree <- (pretty(bombus.par$rare.degree, n=10))
axis.degree <-  standardize.axis(labs.degree,
                                  bombus.par$rare.degree)


## sp.par <- spec.orig[spec.orig$WeightsPar==1, ]
## ## rare.degree (logged)
## labs.degree <- (pretty(sp.par$rare.degree, n=10))
## axis.degree <-  standardize.axis(labs.degree,
##                                   sp.par$rare.degree)


## ***********************************************************************
## bee community diversity and abundance and parasitism
## ***********************************************************************
load(file="saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat_bombus_all.Rdata")
fit.bombus <- fit.parasite

load(file="saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat_social_species.Rdata")
fit.bombus.apicystis <- fit.parasite

load(file="saved/parasiteFit_apis_CrithidiaPresenceApicystisSpp_lat_apis_abundance.Rdata")
fit.apis <- fit.parasite
## https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/
## Generate newdata draws

bombus.cond.effects <- conditional_effects(fit.bombus)

apicystis.cond.effects <- conditional_effects(fit.bombus.apicystis)

apis.cond.effects <- conditional_effects(fit.apis)


## ***************************************************************************
## Crithidia ~ bee diversity

crithidia_beediv <-
  bombus.cond.effects[["CrithidiaPresence.CrithidiaPresence_Net_BeeDiversity"]]

## Add a column to identify from which model these outputs come from.
crithidia_beediv <- mutate(crithidia_beediv, Bee = "Bombus")

# crithidia_beediv_apis <-
#   apis.cond.effects[["CrithidiaPresence.CrithidiaPresence_Net_BeeDiversity"]]
# 
# ## The apis model doesn't have the GenusSpecies random effect thus you need to 
# ## add a column for it to match with bombus.
# crithidia_beediv_apis <- mutate(crithidia_beediv_apis, Bee = "Apis", GenusSpecies = NA)
# 
# ## Join both bombus and apis model outputs in the same dataframe.
# crithidia_beediv<- rbind(crithidia_beediv, crithidia_beediv_apis)

p1.parasite <- ggplot(crithidia_beediv, aes(x = Net_BeeDiversity, y= estimate__)) +
  geom_line(aes(x = Net_BeeDiversity, y= estimate__ , color = Bee), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Bee), alpha=0.4)+
  scale_fill_manual(values =  "#3182bd", 
                    labels =  "Bombus 0.95") +
  scale_color_manual(values =  "#3182bd") +
  labs(x = "Bee diversity", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  guides(color = F, fill = F)+
  scale_x_continuous(
    breaks = axis.bee.div,
    labels =  labs.bee.div) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        text = element_text(size=16))+
  geom_jitter(data=spec.uni,
              aes(y= CrithidiaParasitismRate, x=Net_BeeDiversity),
              width=0.05) 
    
ggsave(p1.parasite, file="figures/parasite_beediv_Crithidia.pdf",
           height=5, width=10)

################################################################################
## Apicystis ~ bee diversity Bombus
################################################################################

apicystis_beediv <-apicystis.cond.effects[["ApicystisSpp.ApicystisSpp_Net_BeeDiversity"]]

#apicystis_beediv <- mutate(apicystis_beediv, Bee = "Bombus")

p2.parasite <- ggplot(apicystis_beediv, aes(x = Net_BeeDiversity, y = estimate__)) +
  geom_line(aes(x = Net_BeeDiversity, y= estimate__), linewidth = 1.5, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__),fill= "#3182bd", alpha=0.4)+
  scale_fill_manual(labels ="Bombus 0.95") +
  labs(x = "Bee diversity", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.bee.div,
    labels =  labs.bee.div) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16))+
  geom_jitter(data=spec.uni,
              aes(y= ApicystisParasitismRate, x=Net_BeeDiversity,
              ), width=0.05)

ggsave(p2.parasite, file="figures/parasite_beediv_Apicystis.pdf",
       height=5, width=10)

## ***************************************************************************
## Crithidia ~ floral diversity

crithidia_floraldiv <-bombus.cond.effects[["CrithidiaPresence.CrithidiaPresence_MeanFloralDiversity"]]
crithidia_floraldiv <- mutate(crithidia_floraldiv, Bee = "Bombus")

# crithidia_floraldiv_apis <-
#   apis.cond.effects[["CrithidiaPresence.CrithidiaPresence_MeanFloralDiversity"]]
# 
# ## The apis model doesn't have the GenusSpecies random effect thus you need to 
# ## add a column for it to match with bombus.
# crithidia_floraldiv_apis <- mutate(crithidia_floraldiv_apis, Bee = "Apis", GenusSpecies = NA)
# 
# ## Join both bombus and apis model outputs in the same dataframe.
# crithidia_floraldiv<- rbind(crithidia_floraldiv, crithidia_floraldiv_apis)

p3.parasite <- ggplot(crithidia_floraldiv, aes(x = MeanFloralDiversity, y= estimate__)) +
  geom_line(aes(x = MeanFloralDiversity, y= estimate__ ,
                color = Bee), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Bee), alpha=0.4) +
  scale_fill_manual(values =  "#3182bd", 
                    labels = "Bombus 0.95") +
  scale_color_manual(values = "#3182bd") +
    labs(x = "Floral diversity", y = "Crithidia prevalence",
         fill = "Credible interval") +
    theme_ms() +
    #theme(legend.position = "bottom") +
  guides(color = FALSE, fill = FALSE)+
    scale_x_continuous(
        breaks = axis.flower.div,
        labels =  labs.flower.div) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
  geom_jitter(data=spec.uni,
              aes(y= ApicystisParasitismRate, x=MeanFloralDiversity,
              ), width=0.05)

ggsave(p3.parasite, file="figures/parasite_floraldiv_Crithidia.pdf",
       height=5, width=10)

## ***************************************************************************
## Apicystis ~ floral diversity

apicystis_floraldiv <- 
  apicystis.cond.effects[["ApicystisSpp.ApicystisSpp_MeanFloralDiversity"]]
#apicystis_floraldiv <- mutate(apicystis_floraldiv, Bee = "Bombus")

p4.parasite <- ggplot(apicystis_floraldiv, aes(x = MeanFloralDiversity, y = estimate__)) +
  geom_line(aes(x = MeanFloralDiversity, y= estimate__), size = 1.5, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha=0.4, fill = "#3182bd")+
  scale_fill_manual(labels ="Bombus 0.95")+
  labs(x = "Floral diversity", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.flower.div,
    labels =  labs.flower.div) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16))+
  geom_jitter(data=spec.uni,
              aes(y= ApicystisParasitismRate, x=MeanFloralDiversity,
              ), width=0.05)

ggsave(p4.parasite, file="figures/parasite_floraldiv_Apicystis.pdf",
       height=5, width=10)
    
    
    
parasite.dilution <- ggarrange(p3.parasite, p1.parasite, p4.parasite, p2.parasite,
                            labels = c("A", "B", "C","D"), 
                            ncol = 2, nrow = 2)

ggsave(parasite.dilution, file="figures/parasite_diversity.pdf", height=8, width=12)

## ***************************************************************************
## crithidia ~ bee abundance

crithidia_beeabun <-
  bombus.cond.effects[["CrithidiaPresence.CrithidiaPresence_Net_BeeAbundance"]]
crithidia_beeabun <- mutate(crithidia_beeabun, Bee = "Bombus")

p5.parasite <- ggplot(crithidia_beeabun, aes(x = Net_BeeAbundance, y = estimate__)) +
  geom_line(aes(x = Net_BeeAbundance, y= estimate__), size = 1.5, ) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha=0.4) +
  #scale_fill_manual(values = "grey", labels ="Bombus 0.95")+
  labs(x = "Bee abundance (log)", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.bee.abund,
    labels =  labs.bee.abund) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
   geom_jitter(data=spec.uni,
               aes(y= CrithidiaParasitismRate, x=Net_BeeAbundance), 
               width=0.05) 


ggsave(p5.parasite, file="figures/crithidia_beeabundance_bombus.pdf",
       height=4, width=6)


################################################################################
## Apicysits ~ bee abundance
################################################################################
apicystis_bombusabun <-
  apicystis.cond.effects[["ApicystisSpp.ApicystisSpp_Net_BombusAbundance"]]
apicystis_bombusabun <- mutate(apicystis_bombusabun, Bee = "Bombus")

p6.parasite <- ggplot(apicystis_bombusabun, aes(x = Net_BombusAbundance, 
                                                y = estimate__)) +
  geom_line(aes(x = Net_BombusAbundance, y= estimate__), size = 1.5, 
            color = "darkgoldenrod3") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Bee), alpha=0.4)+
  scale_fill_manual(values = "darkgoldenrod3", labels ="Bombus 0.95") +
  labs(x = "Bombus abundance (log)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.bombus.abund,
    labels =  labs.bombus.abund) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) + 
  geom_jitter(data=spec.uni,
              aes(y= ApicystisParasitismRate, x=Net_BombusAbundance),width=0.05) 

ggsave(p6.parasite, file="figures/Apicystis_bombusabun_bombus.pdf",
       height=4, width=6)

apicystis_hbabun <-
  apicystis.cond.effects[["ApicystisSpp.ApicystisSpp_Net_HBAbundance"]]
apicystis_hbabun <- mutate(apicystis_hbabun, Bee = "Bombus")

p7.parasite <- ggplot(apicystis_hbabun, aes(x = Net_HBAbundance, 
                                            y = estimate__)) +
  geom_line(aes(x = Net_HBAbundance, y= estimate__), size = 1.5, 
            color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Bee), alpha=0.4)+
  scale_fill_manual(values = "#3182bd", labels ="Bombus 0.95") +
  labs(x = "Apis abundance (log)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.HB.abund,
    labels =  labs.HB.abund) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) + 
  geom_jitter(data=spec.uni,
              aes(y= ApicystisParasitismRate, x=Net_HBAbundance),width=0.05) 

ggsave(p7.parasite, file="figures/Apicystis_hbabun_bombus.pdf",
       height=4, width=6)


parasite.amplification <- ggarrange(p5.parasite, p6.parasite, p7.parasite,
                                    nrow= 1, ncol = 3, heights = 1, widths = c(2, 2, 2),
                                    labels = c("A", "B", "C"))

ggsave(parasite.amplification, file="figures/parasite_amplification.pdf",
       height=4, width=10)

## ***************************************************************************
## parasites ~ diet breadth 
## ***************************************************************************

## ***************************************************************************
## crithidia ~ degree
crithidia_degree <-
  bombus.cond.effects[["CrithidiaPresence.CrithidiaPresence_rare.degree"]]
crithidia_degree <- mutate(crithidia_degree, Bee = "Bombus")

p8.parasite <- ggplot(crithidia_degree, aes(x = rare.degree, y = estimate__)) +
  geom_line(aes(x = rare.degree, y= estimate__), size = 1.5, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Bee), alpha=0.4) +
  scale_fill_manual(values = "#3182bd", labels ="Bombus 0.95")+
    labs(x = "Degree", y = "Crithidia prevalence",
         fill = "Credible interval") +
    theme_ms() +
    #theme(legend.position = "bottom") +
    scale_x_continuous(
        breaks = axis.degree,
        labels =  labs.degree) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
  geom_point(data= bombus.par,
              aes(y= SpCrithidiaParasitismRate, x=rare.degree))
  


ggsave(p8.parasite, file="figures/crithidia_degree.pdf",
       height=5, width=10)

## Apicystis ~ degree
apicystis_degree <-
  apicystis.cond.effects[["ApicystisSpp.ApicystisSpp_rare.degree"]]
apicystis_degree <- mutate(apicystis_degree, Bee = "Bombus")

p9.parasite <- ggplot(apicystis_degree, aes(x = rare.degree, y = estimate__)) +
  geom_line(aes(x = rare.degree, y= estimate__), size = 1.5, color = "grey") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Bee), alpha=0.4) +
  scale_fill_manual(values = "grey", labels ="Bombus 0.95")+
  labs(x = "Degree", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.degree,
    labels =  labs.degree) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
   geom_point(data=bombus.par,
               aes(y= SpApicystisParasitismRate, x=rare.degree)) 

ggsave(p9.parasite, file="figures/apicystis_degree.pdf",
       height=5, width=10)

parasite.traits <- ggarrange(p8.parasite, p9.parasite, nrow=1,
                             labels = c("A", "B"))
ggsave(parasite.traits, file="figures/parasite_traits.pdf",
       height= 6, width=10)


