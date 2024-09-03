rm(list=ls())
## setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')
## setwd("C:/Users/na_ma/Dropbox (University of Oregon)/skyIslands")
setwd('~/Dropbox (University of Oregon)/skyislands')
## Script for plotting all of the important explanatory variables.
library(ggpubr)
library(tidyverse)
library(tidybayes)
setwd("analysis/parasites")
load(file="saved/spec_weights.Rdata")
source("src/misc.R")
source("src/ggplotThemes.R")

## ***********************************************************************
## plotting, unscaling labs
## ***********************************************************************
## lat

spec.uni <- spec.orig[spec.orig$Weights ==1,]
labs.lat.x <- pretty(c(spec.uni$Lat),
                     n=10)
axis.lat.x <-  standardize.axis(labs.lat.x, spec.uni$Lat)
## bloom abundance
labs.bloom.abund <- (pretty(c(spec.uni$MeanFloralAbundance), n=6))
axis.bloom.abund <-  standardize.axis(labs.bloom.abund,
                                      spec.uni$MeanFloralAbundance)
## flower div
labs.flower.div <- (pretty(spec.uni$MeanFloralDiversity, n=5))
axis.flower.div <-  standardize.axis(labs.flower.div,
                                     spec.uni$MeanFloralDiversity)
## HB abund
labs.HB.abund <- (pretty(c(spec.uni$Net_HBAbundance), n=5))
axis.HB.abund <-  standardize.axis(labs.HB.abund, spec.uni$Net_HBAbundance)
## bombus abund
labs.bombus.abund <- (pretty(c(spec.uni$Net_BombusAbundance), n=5))
axis.bombus.abund <-  standardize.axis(labs.bombus.abund, spec.uni$Net_BombusAbundance)
## non hb non bombus abund
labs.bee.abund <- (pretty(c(spec.uni$Net_NonBombusHBAbundance), n=5))
axis.bee.abund <-  standardize.axis(labs.bee.abund, spec.uni$Net_NonBombusHBAbundance)
## bee diversity
labs.bee.div <- (pretty(c(spec.uni$Net_BeeDiversity), n=5))
axis.bee.div <-  standardize.axis(labs.bee.div,
                                  spec.uni$Net_BeeDiversity)

spec.sp <-  spec.orig[spec.orig$WeightsSp ==1 & spec.orig$Genus == "Bombus",]
## meanITD
labs.itd <- (pretty(spec.sp$MeanITD, n=10))
axis.itd <-  standardize.axis(labs.itd,
                              spec.sp$MeanITD)
## rare.degree
labs.degree <- (pretty(spec.sp$rare.degree, n=10))
axis.degree <-  standardize.axis(labs.itd,
                                 spec.sp$rare.degree)

## ***********************************************************************
## bee community diversity and abundance and parasitism
## ***********************************************************************
load(file="../../../skyIslands_saved/parasite-results/saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat_all.Rdata")

# We want the standarized data for the predictions (spec.data)
spec.bombus <- spec.net[spec.net$Genus == "Bombus",]
#bombus.par <- spec.bombus[spec.bombus$WeightsSp == 1,]
data.site <- spec.net[spec.net$Weights == 1,]
## ***********************************************************************
## bee community diversity and latitude
## ***********************************************************************

## Community level visuals

newdata.lat <- crossing(Lat =
                          seq(min(data.site$Lat),
                              max(data.site$Lat),
                              length.out=10),
                        Net_BeeAbundance = 0,
                        MeanFloralDiversity = 0,
                        Net_BeeDiversity = 0,
                        SRDoyPoly1 = 0,
                        SRDoyPoly2 = 0,
                        Area = 0, 
                        Year = "2012",
                        Site = "SC", 
                        GenusSpecies = "Bombus centralis",
                        Weights=1,
                        WeightsPar=1
)

## predict values based on generated data and model parameters
pred_lat <- fit.parasite %>% 
  epred_draws(newdata = newdata.lat,
              resp = "NetBeeDiversity")

## to see range of predicted values
pred_lat %>%
  group_by(Net_BeeDiversity) %>%
  summarise(mean(.epred))

beediv_lat <- ggplot(pred_lat, aes(x = .epred, y = Net_BeeDiversity)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Latitude", y = "Bee Species Diversity",
       fill = "Credible interval") +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  scale_y_continuous(
    breaks = axis.bee.div,
    labels =  labs.bee.div) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  ##theme_dark_black()+
  geom_point(data=data.site,
             aes(y= Net_BeeDiversity, x=Lat),
             color="grey40", cex=2)

ggsave(beediv_lat, file="figures/Lat_beediv.pdf",
       height=4, width=5)

################################################################################
## Plant community diversity and latitude
################################################################################

newdata.lat <- crossing(Lat =
                          seq(min(data.site$Lat),
                              max(data.site$Lat),
                              length.out=10),
                        Net_BeeAbundance = 0,
                        MeanFloralDiversity = 0,
                        Net_BeeDiversity = 0,
                        SRDoyPoly1 = 0,
                        SRDoyPoly2 = 0,
                        Area = 0, 
                        Year = "2012",
                        Site = "SC", 
                        GenusSpecies = "Bombus centralis",
                        Weights=1
)

## predict values based on generated data and model parameters
pred_lat <- fit.parasite %>% 
  epred_draws(newdata = newdata.lat,
              resp = "MeanFloralDiversity")

## to see range of predicted values
pred_lat %>%
  group_by(MeanFloralDiversity) %>%
  summarise(mean(.epred))

plantdiv_lat <- ggplot(pred_lat, aes(x = .epred, y = MeanFloralDiversity)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Latitude", y = "Flowering Species Diversity",
       fill = "Credible interval") +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  scale_y_continuous(
    breaks = axis.flower.div,
    labels =  labs.flower.div) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  ##theme_dark_black()+
  geom_point(data=data.site,
             aes(y= MeanFloralDiversity, x=Lat),
             color="grey40", cex=2)

ggsave(plantdiv_lat, file="figures/Lat_floraldiv.pdf",
       height=4, width=5)
################################################################################
## Plant diversity and bee diversity
################################################################################
newdata.div <- crossing(MeanFloralDiversity =
                          seq(min(data.site$MeanFloralDiversity),
                              max(data.site$MeanFloralDiversity),
                              length.out=10),
                        Net_BeeAbundance = 0,
                        Net_BeeDiversity = 0,
                        SRDoyPoly1 = 0,
                        SRDoyPoly2 = 0,
                        Area = 0, 
                        Year = "2012",
                        Site = "SC", 
                        GenusSpecies = "Bombus centralis",
                        Weights=1
)

## predict values based on generated data and model parameters
pred_div <- fit.parasite %>% 
  epred_draws(newdata = newdata.div,
              resp = "NetBeeDiversity")

## to see range of predicted values
pred_div %>%
  group_by(NetBeeDiversity) %>%
  summarise(mean(.epred))

plantdiv_beediv <- ggplot(pred_div, aes(x = .epred, y = MeanFloralDiversity)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Bee Species Diversity", y = "Floral Diversity",
       fill = "Credible interval") +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.bee.div,
    labels =  labs.bee.div) +
  scale_y_continuous(
    breaks = axis.flower.div,
    labels =  labs.flower.div) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  ##theme_dark_black()+
  geom_point(data=data.site,
             aes(y= Net_BeeDiversity, x= MeanFloralDiversity),
             color="grey40", cex=2)

ggsave(plantdiv_beediv, file="figures/Beediv_floraldiv.pdf",
       height=4, width=5)

################################################################################
## Bee abundance and lat
################################################################################
newdata.lat <- crossing(Lat =
                          seq(min(data.site$Lat),
                              max(data.site$Lat),
                              length.out=10),
                        Net_BeeAbundance = 0,
                        MeanFloralDiversity = 0,
                        Net_BeeDiversity = 0,
                        SRDoyPoly1 = 0,
                        SRDoyPoly2 = 0,
                        Area = 0, 
                        Year = "2012",
                        Site = "SC", 
                        GenusSpecies = "Bombus centralis",
                        Weights=1
)

## predict values based on generated data and model parameters
pred_lat <- fit.parasite %>% 
  epred_draws(newdata = newdata.lat,
              resp = "NetBombusAbundance")

## to see range of predicted values
pred_lat %>%
  group_by(Net_BombusAbundance) %>%
  summarise(mean(.epred))

bombusabun_lat <- ggplot(pred_lat, aes(x = .epred, y = Net_BombusAbundance)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Latitude", y = "Bombus Abudance",
       fill = "Credible interval") +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  scale_y_continuous(
    breaks = axis.bombus.abund,
    labels =  labs.bombus.abund) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  ##theme_dark_black()+
  geom_point(data=data.site,
             aes(y= Net_BombusAbundance, x=Lat),
             color="grey40", cex=2)

ggsave(bombusabun_lat, file="figures/Lat_bombus_abudance.pdf",
       height=4, width=5)

## Honeybee abundance
newdata.lat <- crossing(Lat =
                          seq(min(data.site$Lat),
                              max(data.site$Lat),
                              length.out=10),
                        Net_BeeAbundance = 0,
                        MeanFloralDiversity = 0,
                        Net_BeeDiversity = 0,
                        SRDoyPoly1 = 0,
                        SRDoyPoly2 = 0,
                        Area = 0, 
                        Year = "2012",
                        Site = "SC", 
                        GenusSpecies = "Bombus centralis",
                        Weights=1
)

## predict values based on generated data and model parameters
pred_lat <- fit.parasite %>% 
  epred_draws(newdata = newdata.lat,
              resp = "NetHBAbundance")

## to see range of predicted values
pred_lat %>%
  group_by(Net_BombusAbundance) %>%
  summarise(mean(.epred))

HBabun_lat <- ggplot(pred_lat, aes(x = .epred, y = Net_HBAbundance)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Latitude", y = "Bombus Abudance",
       fill = "Credible interval") +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  scale_y_continuous(
    breaks = axis.HB.abund,
    labels =  labs.HB.abund) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  ##theme_dark_black()+
  geom_point(data=data.site,
             aes(y= Net_HBAbundance, x=Lat),
             color="grey40", cex=2)

ggsave(HBabun_lat, file="figures/Lat_HB_abudance.pdf",
       height=4, width=5)

################################################################################
##SRDoyPoly1 
################################################################################
newdata.srdoypol1 <- crossing(SRDoyPoly1 =
                          seq(min(data.site$SRDoyPoly1),
                              max(data.site$SRDoyPoly1),
                              length.out=10),
                        Net_BeeAbundance = 0,
                        MeanFloralDiversity = 0,
                        Net_BeeDiversity = 0,
                        Lat = 0,
                        SRDoyPoly2 = 0,
                        Area = 0, 
                        Year = "2012",
                        Site = "SC", 
                        GenusSpecies = "Bombus centralis",
                        Weights=1
)

## predict values based on generated data and model parameters
pred_SRDoyPoly1 <- fit.parasite %>% 
  epred_draws(newdata = newdata.SRDoyPoly1,
              resp = "NetBeeDiversity")

## to see range of predicted values
pred_SRDoyPoly1 %>%
  group_by(Net_BeeDiversity) %>%
  summarise(mean(.epred))

beediv_lat <- ggplot(pred_SRDoyPoly1, aes(x = .epred, y = Net_BeeDiversity)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "First Polynomial Term", y = "Bee Species Diversity",
       fill = "Credible interval") +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.,
    labels =  labs.lat.x) +
  scale_y_continuous(
    breaks = axis.bee.div,
    labels =  labs.bee.div) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  ##theme_dark_black()+
  geom_point(data=data.site,
             aes(y= Net_BeeDiversity, x=Lat),
             color="grey40", cex=2)

ggsave(beediv_lat, file="figures/Lat_beediv.pdf",
       height=4, width=5)

################################################################################
## Plant community diversity and latitude
################################################################################

newdata.lat <- crossing(Lat =
                          seq(min(data.site$Lat),
                              max(data.site$Lat),
                              length.out=10),
                        Net_BeeAbundance = 0,
                        MeanFloralDiversity = 0,
                        Net_BeeDiversity = 0,
                        SRDoyPoly1 = 0,
                        SRDoyPoly2 = 0,
                        Area = 0, 
                        Year = "2012",
                        Site = "SC", 
                        GenusSpecies = "Bombus centralis",
                        Weights=1
)

## predict values based on generated data and model parameters
pred_lat <- fit.parasite %>% 
  epred_draws(newdata = newdata.lat,
              resp = "MeanFloralDiversity")

## to see range of predicted values
pred_lat %>%
  group_by(MeanFloralDiversity) %>%
  summarise(mean(.epred))

plantdiv_lat <- ggplot(pred_lat, aes(x = .epred, y = MeanFloralDiversity)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Latitude", y = "Flowering Species Diversity",
       fill = "Credible interval") +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  scale_y_continuous(
    breaks = axis.flower.div,
    labels =  labs.flower.div) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  ##theme_dark_black()+
  geom_point(data=data.site,
             aes(y= MeanFloralDiversity, x=Lat),
             color="grey40", cex=2)

ggsave(plantdiv_lat, file="figures/Lat_floraldiv.pdf",
       height=4, width=5)
################################################################################
## Plant diversity and bee diversity
################################################################################
newdata.div <- crossing(MeanFloralDiversity =
                          seq(min(data.site$MeanFloralDiversity),
                              max(data.site$MeanFloralDiversity),
                              length.out=10),
                        Net_BeeAbundance = 0,
                        Net_BeeDiversity = 0,
                        SRDoyPoly1 = 0,
                        SRDoyPoly2 = 0,
                        Area = 0, 
                        Year = "2012",
                        Site = "SC", 
                        GenusSpecies = "Bombus centralis",
                        Weights=1
)

## predict values based on generated data and model parameters
pred_div <- fit.parasite %>% 
  epred_draws(newdata = newdata.div,
              resp = "NetBeeDiversity")

## to see range of predicted values
pred_div %>%
  group_by(NetBeeDiversity) %>%
  summarise(mean(.epred))

plantdiv_beediv <- ggplot(pred_div, aes(x = .epred, y = MeanFloralDiversity)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Bee Species Diversity", y = "Floral Diversity",
       fill = "Credible interval") +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.bee.div,
    labels =  labs.bee.div) +
  scale_y_continuous(
    breaks = axis.flower.div,
    labels =  labs.flower.div) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  ##theme_dark_black()+
  geom_point(data=data.site,
             aes(y= Net_BeeDiversity, x= MeanFloralDiversity),
             color="grey40", cex=2)

ggsave(plantdiv_beediv, file="figures/Beediv_floraldiv.pdf",
       height=4, width=5)

################################################################################
## Bee abundance and lat
################################################################################
newdata.lat <- crossing(Lat =
                          seq(min(data.site$Lat),
                              max(data.site$Lat),
                              length.out=10),
                        Net_BeeAbundance = 0,
                        MeanFloralDiversity = 0,
                        Net_BeeDiversity = 0,
                        SRDoyPoly1 = 0,
                        SRDoyPoly2 = 0,
                        Area = 0, 
                        Year = "2012",
                        Site = "SC", 
                        GenusSpecies = "Bombus centralis",
                        Weights=1
)

## predict values based on generated data and model parameters
pred_lat <- fit.parasite %>% 
  epred_draws(newdata = newdata.lat,
              resp = "NetBombusAbundance")

## to see range of predicted values
pred_lat %>%
  group_by(Net_BombusAbundance) %>%
  summarise(mean(.epred))

bombusabun_lat <- ggplot(pred_lat, aes(x = .epred, y = Net_BombusAbundance)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Latitude", y = "Bombus Abudance",
       fill = "Credible interval") +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  scale_y_continuous(
    breaks = axis.bombus.abund,
    labels =  labs.bombus.abund) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  ##theme_dark_black()+
  geom_point(data=data.site,
             aes(y= Net_BombusAbundance, x=Lat),
             color="grey40", cex=2)

ggsave(bombusabun_lat, file="figures/Lat_bombus_abudance.pdf",
       height=4, width=5)
