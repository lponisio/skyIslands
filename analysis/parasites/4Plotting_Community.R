rm(list=ls())
## setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')
## setwd("C:/Users/na_ma/Dropbox (University of Oregon)/skyIslands")
setwd('~/Dropbox (University of Oregon)/skyislands')
## Script for plotting all of the important explanatory variables.
library(ggpubr)
library(tidyverse)
library(tidybayes)
library(patchwork)
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


## ***********************************************************************
## bee community diversity and abundance and parasitism
## ***********************************************************************
load(file="../../../skyIslands_saved/parasite-results/saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat_all.Rdata")
fit.bombus <- fit.parasite

cond.effects <- conditional_effects(fit.bombus)
## Community level visuals
## ***********************************************************************
## bee community diversity and latitude
## ***********************************************************************

lat_beediv <-
  cond.effects[["NetBeeDiversity.NetBeeDiversity_Lat"]]

p1 <- ggplot(lat_beediv, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha=0.4, fill = "#e6550d") +
  labs(x = "Latitude", y = "Bee diversity",
       fill = "Credible interval")+
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  scale_y_continuous(
    breaks = axis.bee.div,
    labels =  labs.bee.div) +
  theme(axis.title.x = element_blank(),#Remember to change the x axis label when using this graph only
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  geom_point(data= spec.uni,
             aes(y= Net_BeeDiversity, x=Lat),
             color="grey40", cex=2)

ggsave(p1, file="figures/Lat_beediv.pdf",
       height=4, width=5)

################################################################################
## Plant community diversity and latitude
################################################################################

lat_floraldiv <-
  cond.effects[["MeanFloralDiversity.MeanFloralDiversity_Lat"]]

p2 <- ggplot(lat_floraldiv, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha=0.4, fill = "darkgoldenrod3") +
  labs(x = "Latitude", y = "Mean floral diversity",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  scale_y_continuous(
    breaks = axis.flower.div,
    labels =  labs.flower.div) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  ##theme_dark_black()+
  geom_point(data=spec.uni,
             aes(y= MeanFloralDiversity, x=Lat),
             color="grey40", cex=2)

ggsave(p2, file="figures/Lat_floraldiv.pdf",
       height=4, width=5)


################################################################################
## Bee abundance and lat
################################################################################

lat_bombusabund <-
  cond.effects[["NetBombusAbundance.NetBombusAbundance_Lat"]]


p3 <- ggplot(lat_bombusabund, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha=0.4, fill = "#e6550d") +
  labs(x = "Latitude (log)", y = "Bombus abundance (log)",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  scale_y_continuous(
    breaks = axis.bombus.abund,
    labels =  labs.bombus.abund) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  ##theme_dark_black()+
  geom_point(data=spec.uni,
             aes(y= Net_BombusAbundance, x=Lat),
             color="grey40", cex=2)

ggsave(p3, file="figures/Lat_bombus_abudance.pdf",
       height=4, width=5)

## Honeybee abundance
lat_apisabund <-
  cond.effects[["NetHBAbundance.NetHBAbundance_Lat"]]

p4 <- ggplot(lat_apisabund, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha=0.4, fill = "#e6550d") +
  labs(x = "Latitude (log)", y = "Apis abundance (log)",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  scale_y_continuous(
    breaks = axis.HB.abund,
    labels =  labs.HB.abund) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  ##theme_dark_black()+
  geom_point(data=spec.uni,
             aes(y= Net_HBAbundance, x=Lat),
             color="grey40", cex=2)

ggsave(p4, file="figures/Lat_HB_abudance.pdf",
       height=4, width=5)

lat_community<- ggarrange(p1,p2,p3,p4, #plots that are going to be included in this multipanel figure
                       labels = c("A", "B", "C","D"), #labels given each panel 
                       ncol = 2, nrow = 2 #adjust plot space 
                       )
ggsave(lat_community, file="figures/lat_community.pdf",
       height=8, width=12)
################################################################################
## Plant diversity and bee diversity
################################################################################
beediv_floraldiv <-
  cond.effects[["NetBeeDiversity.NetBeeDiversity_MeanFloralDiversity"]]

plantdiv_beediv <- ggplot(beediv_floraldiv, aes(x = MeanFloralDiversity, y = estimate__)) +
  geom_line(aes(x = MeanFloralDiversity, y= estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha=0.4, fill = "#3182bd") +
  labs(y = "Bee Species Diversity", x = "Floral Diversity",
       fill = "Credible interval") +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.flower.div,
    labels =  labs.flower.div) +
  scale_y_continuous(
    breaks = axis.bee.div,
    labels =  labs.bee.div) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  ##theme_dark_black()+
  geom_point(data=spec.uni,
             aes(y= Net_BeeDiversity, x= MeanFloralDiversity),
             color="grey40", cex=2)

ggsave(plantdiv_beediv, file="figures/Beediv_floraldiv.pdf",
       height=4, width=5)

