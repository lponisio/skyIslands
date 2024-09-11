rm(list=ls())
## setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')
 setwd("C:/Users/na_ma/Dropbox (University of Oregon)/skyIslands")
##setwd('~/Dropbox (University of Oregon)/skyislands')
## Script for plotting all of the important explanatory variables.
library(ggpubr)
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
load(file="../../../skyIslands_saved/parasite-results/saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat_all.Rdata")
fit.bombus <- fit.parasite

(file="../../../skyIslands_saved/parasite-results/saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat_ba.Rdata")
fit.bombus.ba <- fit.parasite

(file="../../../skyIslands_saved/parasite-results/saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat_ss.Rdata")
fit.bombus.ss <- fit.parasite

load(file="../../../skyIslands_saved/parasite-results/saved/parasiteFit_apis_CrithidiaPresenceApicystisSpp_lat_ss.Rdata")
fit.apis <- fit.parasite
## https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/
## Generate newdata draws

bombus.cond.effects <- conditional_effects(fit.bombus)
bombus.cond.effects2<- conditional_effects(fit.bombus.ss)
bombus.cond.effects3<- conditional_effects(fit.bombus.ba)
apis.cond.effects <- conditional_effects(fit.apis)


## ***************************************************************************
## Crithidia ~ bee diversity

crithidia_beediv <-
  bombus.cond.effects[["CrithidiaPresence.CrithidiaPresence_Net_BeeDiversity"]]

## Add a column to identify from which model these outputs come from.
crithidia_beediv <- mutate(crithidia_beediv, Bee = "Bombus")

crithidia_beediv_apis <-
  apis.cond.effects[["CrithidiaPresence.CrithidiaPresence_Net_BeeDiversity"]]

## The apis model doesn't have the GenusSpecies random effect thus you need to 
## add a column for it to match with bombus.
crithidia_beediv_apis <- mutate(crithidia_beediv_apis, Bee = "Apis", GenusSpecies = NA)

## Join both bombus and apis model outputs in the same dataframe.
crithidia_beediv<- rbind(crithidia_beediv, crithidia_beediv_apis)

p1.parasite <- ggplot(crithidia_beediv, aes(x = Net_BeeDiversity, y= estimate__)) +
  geom_line(aes(x = Net_BeeDiversity, y= estimate__ , color = Bee), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Bee), alpha=0.4)+
  scale_fill_manual(values = c("darkgoldenrod3", "#3182bd"), 
                    labels = c("Apis 0.95", "Bombus 0.95")) +
  scale_color_manual(values = c("darkgoldenrod3", "#3182bd")) +
  labs(x = "Bee diversity", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.bee.div,
    labels =  labs.bee.div) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        text = element_text(size=16))+
     geom_jitter(data=spec.uni,
             aes(y= CrithidiaParasitismRate, x=Net_BeeDiversity,
                 ), width=0.05) 
    
ggsave(p1.parasite, file="figures/parasite_beediv_Crithidia.pdf",
           height=5, width=10)

################################################################################
## Apicystis ~ bee diversity Bombus
################################################################################
apicystis_beediv <-bombus.cond.effects[["ApicystisSpp.ApicystisSpp_Net_BeeDiversity"]]

apicystis_beediv <- mutate(apicystis_beediv, Bee = "Bombus")

p2.parasite <- ggplot(apicystis_beediv, aes(x = Net_BeeDiversity, y = estimate__)) +
  geom_line(aes(x = Net_BeeDiversity, y= estimate__), size = 1.5, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill= Bee), alpha=0.4)+
  scale_fill_manual(values = "#3182bd", labels ="Bombus 0.95") +
  labs(x = "Bee diversity", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "bottom") +
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

p3.parasite <- ggplot(crithidia_floraldiv, aes(x = MeanFloralDiversity, y = estimate__)) +
  geom_line(aes(x = MeanFloralDiversity, y= estimate__), size = 1.5, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Bee), alpha=0.4) +
  scale_fill_manual(values = "#3182bd", labels ="Bombus 0.95") +
    labs(x = "Floral diversity", y = "Crithidia prevalence",
         fill = "Credible interval") +
    theme_ms() +
    theme(legend.position = "bottom") +
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
  bombus.cond.effects2[["ApicystisSpp.ApicystisSpp_MeanFloralDiversity"]]
apicystis_floraldiv <- mutate(apicystis_floraldiv, Bee = "Bombus")

p4.parasite <- ggplot(crithidia_floraldiv, aes(x = MeanFloralDiversity, y = estimate__)) +
  geom_line(aes(x = MeanFloralDiversity, y= estimate__), size = 1.5, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Bee), alpha=0.4)+
  scale_fill_manual(values = "#3182bd", labels ="Bombus 0.95")+
  labs(x = "Floral diversity", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "bottom") +
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
    
    
    
parasite.dilution <- ggarrange(p1.parasite, p3.parasite,p2.parasite, p4.parasite, 
                            labels = c("A", "B", "C","D"), 
                            ncol = 2, nrow = 2,
                          common.legend = T,
                          legend = "bottom")
ggsave(parasite.dilution, file="figures/parasite_diversity.pdf", height=8, width=12)

## ***************************************************************************
## crithidia ~ bee abundance

crithidia_beeabun <-
  bombus.cond.effects[["CrithidiaPresence.CrithidiaPresence_Net_BeeAbundance"]]
crithidia_beeabun <- mutate(crithidia_beeabun, Bee = "Bombus")

p5.parasite <- ggplot(crithidia_beeabun, aes(x = Net_BeeAbundance, y = estimate__)) +
  geom_line(aes(x = Net_BeeAbundance, y= estimate__), size = 1.5, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Bee), alpha=0.4) +
  scale_fill_manual(values = "#3182bd", labels ="Bombus 0.95")+
  labs(x = "Bee abundance (log)", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "bottom") +
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

## crithidia ~ bombus abundance 
## Note this is using the bombus abundance model only

 crithidia_bombusabun <-
   bombus.cond.effects3[["CrithidiaPresence.CrithidiaPresence_Net_BombusAbundance"]]
# crithidia_bombusabun <- mutate(crithidia_beeabun, Bee = "Bombus")
 
 p6.parasite <- ggplot(crithidia_bombusabun, aes(x = Net_BombusAbundance, y = estimate__)) +
   geom_line(aes(x = Net_BombusAbundance, y= estimate__), size = 1.5, color = "#3182bd") +
   geom_ribbon(aes(ymin = lower__, ymax = upper__),fill ="#3182bd", alpha=0.4) +
   scale_fill_manual(labels ="Bombus 0.95")+
   labs(x = "Bombus abundance (log)", y = "Crithidia prevalence",
        fill = "Credible interval") +
   theme_ms() +
theme(legend.position = "bottom") +
   scale_x_continuous(
     breaks = axis.bombus.abund,
     labels =  labs.bombus.abund) +
   theme(axis.title.x = element_text(size=16),
         axis.title.y = element_text(size=16),
text = element_text(size=16)) +
   geom_jitter(data=spec.uni,
               aes(y= CrithidiaParasitismRate, x=Net_BombusAbundance), 
               width=0.05) 
 
 
 ggsave(p6.parasite, file="figures/crithidia_bombusabundance_bombus.pdf",
        height=4, width=6)
################################################################################
## Apicysits ~ bombus abundance
################################################################################
apicystis_beeabun <-
  bombus.cond.effects2[["ApicystisSpp.ApicystisSpp_Net_BombusAbundance"]]
apicystis_beeabun <- mutate(apicystis_beeabun, Bee = "Bombus")

p7.parasite <- ggplot(apicystis_beeabun, aes(x = Net_BombusAbundance, y = estimate__)) +
  geom_line(aes(x = Net_BombusAbundance, y= estimate__), size = 1, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Bee), alpha=0.4)+
  scale_fill_manual(values = "#3182bd", labels ="Bombus 0.95") +
  labs(x = "Bombus abundance (log)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.bombus.abund,
    labels =  labs.bombus.abund) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) + 
  geom_jitter(data=spec.uni,
              aes(y= ApicystisParasitismRate, x=Net_BombusAbundance),width=0.05) 

ggsave(p7.parasite, file="figures/Apicystis_bombusabun_bombus.pdf",
       height=4, width=6)

## Apicystis ~ HB abundance
apicystis_HBabun <-
  bombus.cond.effects2[["ApicystisSpp.ApicystisSpp_Net_HBAbundance"]]
#apicystis_HBabun <- mutate(apicystis_beeabun, Bee = "Bombus")

p8.parasite <- ggplot(apicystis_HBabun, aes(x = Net_HBAbundance, y = estimate__)) +
  geom_line(aes(x = Net_HBAbundance, y= estimate__), size = 1, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill="#3182bd", alpha=0.4)+
  scale_fill_manual(labels ="Bombus 0.95") +
  labs(x = "Apis abundance (log)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.HB.abund,
    labels =  labs.HB.abund) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) + 
  geom_jitter(data=spec.uni,
              aes(y= ApicystisParasitismRate, x=Net_HBAbundance),width=0.05) 

ggsave(p8.parasite, file="figures/Apicystis_HBabun_bombus.pdf",
       height=4, width=6)

parasite.amplification <- ggarrange(p5.parasite, p6.parasite, p7.parasite, p8.parasite,
                                    nrow= 2, ncol = 2,
                                    labels = c("A", "B", "C", "D"),
                                    common.legend = TRUE,
                                    legend="bottom")

ggsave(parasite.amplification, file="figures/parasite_amplification.pdf",
       height=6, width=10)

## ***************************************************************************
## parasites ~ diet breadth 
## ***************************************************************************

## ***************************************************************************
## crithidia ~ degree
crithidia_degree <-
  bombus.cond.effects[["CrithidiaPresence.CrithidiaPresence_rare.degree"]]
crithidia_degree <- mutate(crithidia_degree, Bee = "Bombus")

p9.parasite <- ggplot(crithidia_degree, aes(x = rare.degree, y = estimate__)) +
  geom_line(aes(x = rare.degree, y= estimate__), size = 1.5, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Bee), alpha=0.4) +
  scale_fill_manual(values = "#3182bd", labels ="Bombus 0.95")+
    labs(x = "Degree", y = "Crithidia prevalence",
         fill = "Credible interval") +
    theme_ms() +
    theme(legend.position = "bottom") +
    scale_x_continuous(
        breaks = axis.degree,
        labels =  labs.degree) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
  geom_point(data= bombus.par,
              aes(y= SpCrithidiaParasitismRate, x=rare.degree))
  


ggsave(p9.parasite, file="figures/crithidia_degree.pdf",
       height=5, width=10)

## Apicystis ~ degree
apicystis_degree <-
  bombus.cond.effects2[["ApicystisSpp.ApicystisSpp_rare.degree"]]
apicystis_degree <- mutate(apicystis_degree, Bee = "Bombus")

p10.parasite <- ggplot(apicystis_degree, aes(x = rare.degree, y = estimate__)) +
  geom_line(aes(x = rare.degree, y= estimate__), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha=0.4) +
  scale_fill_manual(labels ="Bombus 0.95") +
  labs(x = "Degree", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.degree,
    labels =  labs.degree) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
   geom_point(data=bombus.par,
               aes(y= SpApicystisParasitismRate, x=rare.degree)) 

ggsave(p8.parasite, file="figures/apicystis_degree.pdf",
       height=5, width=10)

parasite.traits <- ggarrange(p9.parasite, p10.parasite, nrow=1,
                             labels = c("A", "B"),
                             common.legend = TRUE,
                             legend="bottom")
ggsave(parasite.traits, file="figures/parasite_traits.pdf",
       height= 6, width=10)
## ***********************************************************************
## lat and crithidia
## ***********************************************************************

crithidia_lat <-
  bombus.cond.effects[["CrithidiaPresence.CrithidiaPresence_Lat"]]
crithidia_lat <- mutate(crithidia_lat, Bee = "Bombus")

p11.parasite <- ggplot(crithidia_lat, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size = 1.5, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Bee), alpha=0.4) +
  scale_fill_manual(values = "#3182bd", labels ="Bombus 0.95") +
  labs(x = "Latitude (log)", y = "Crithidia prevalence",
       title = "Bombus") +
  ggtitle("Bombus")+
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        text = element_text(size=16),
        plot.title = element_text(color = "black")) +
  geom_point(data= spec.uni,aes(y= CrithidiaParasitismRate, x= Lat), cex=2)



ggsave(p11.parasite, file="figures/Lat_crithidia.pdf",
       height=4, width=5)

## ***********************************************************************
## lat and apicystis
## ***********************************************************************

apicystis_lat <-
  bombus.cond.effects2[["ApicystisSpp.ApicystisSpp_Lat"]]
apicystis_lat <- mutate(apicystis_lat, Bee = "Bombus")

p12.parasite <- ggplot(apicystis_lat, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size = 1.5, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Bee), alpha=0.4) +
  scale_fill_manual(values = "#3182bd", labels ="Bombus 0.95") +
  labs(x = "Latitude (log)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16))+
   geom_jitter(data= spec.uni,
               aes(y= ApicystisParasitismRate, x= Lat), cex=2) 



ggsave(p12.parasite, file="figures/Lat_apicystis.pdf",
       height=4, width=5)

################################################################################
## crithida ~ lat in Apis
crithidia_lat_apis <-
  apis.cond.effects[["CrithidiaPresence.CrithidiaPresence_Lat"]]
crithidia_lat_apis <- mutate(crithidia_lat_apis, Bee = "Apis")

p13.parasite <- ggplot(crithidia_lat_apis, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size= 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha=0.4) +
  scale_fill_manual(labels ="Apis 0.95") +
  labs(x = "Latitude (log)", y = "Crithidia prevalence",
       title = "Apis") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        text = element_text(size=16), 
        plot.title = element_text(color = "black")) +
geom_jitter(data= spec.uni,
             aes(y= CrithidiaParasitismRate, x= Lat), cex=2) 



ggsave(p13.parasite, file="figures/Lat_crithidia_apis.pdf",
       height=4, width=5)

################################################################################
## Lat and apicystis in apis
################################################################################
apicystis_lat_apis <-
  apis.cond.effects[["ApicystisSpp.ApicystisSpp_Lat"]]
apicystis_lat_apis <- mutate(apicystis_lat_apis, Bee = "Apis")

p14.parasite <- ggplot(apicystis_lat_apis, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha=0.4) +
  scale_fill_manual(labels ="Apis 0.95") +
  labs(x = "Latitude (log)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
geom_jitter(data=spec.uni,
             aes(y= ApicystisParasitismRate, x= Lat), cex=2) 



ggsave(p14.parasite, file="figures/Lat_apicystis_apis.pdf",
       height=4, width=5)



lat_parasite<- ggarrange(p11.parasite, p14.parasite, p12.parasite, p13.parasite,#plots that are going to be included in this multipanel figure
                          labels = c("A", "C", "B","D"), #labels given each panel 
                          ncol = 2, nrow = 2)
 
ggsave(lat_parasite, file="figures/lat_parasites.pdf",
       height=8, width=12)


