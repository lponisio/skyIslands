## setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')
setwd("C:/Users/na_ma/Dropbox (University of Oregon)/skyIslands")
## setwd('~/Dropbox (University of Oregon)/skyislands')
## Script for plotting all of the important explanatory variables.

setwd("analysis/parasites")
rm(list=ls())
source("src/misc.R")
source("src/writeResultsTable.R")
source("src/makeMultiLevelData.R")
source("src/runParasiteModels.R")
source("src/standardize_weights.R")
source("src/ggplotThemes.R")
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
source("src/init.R")

spec.orig <- spec.net
site.orig <- site.sum

## ***********************************************************************
## plotting, unscaling labs
## ***********************************************************************
## lat

spec.orig$Lat <- log(spec.orig$Lat)
labs.lat.x <- (pretty(c(spec.orig$Lat),
                      n=10))
axis.lat.x <-  standardize.axis(labs.lat.x, spec.orig$Lat)
## bloom abundance
labs.bloom.abund <- (pretty(c(0, spec.orig$MeanFloralAbundance), n=6))
axis.bloom.abund <-  standardize.axis(labs.bloom.abund,
                                      spec.orig$MeanFloralAbundance)
## area (notcurrently used due to colinearity with lat)
spec.orig$Area <- log(spec.orig$Area)
labs.area <- (pretty(spec.orig$Area, n=6))
axis.area <-  standardize.axis(labs.area,
                               spec.orig$Area)
## bee abund
## labs.bee.abund2 <- (pretty(c(0, spec.orig$PollAbundance),n=6))
## axis.bee.abund2 <- labs.bee.abund2
## not standardized so can use poisson
## axis.bee.abund2 <-  standardize.axis(labs.bee.abund2,
##                                       spec.orig$PollAbundance)

## flower div
labs.flower.div <- (pretty(spec.orig$MeanFloralDiversity, n=5))
axis.flower.div <-  standardize.axis(labs.flower.div,
                                     spec.orig$MeanFloralDiversity)
## bee abund
## labs.bee.abund <- (pretty(c(0, spec.orig$PollAbundance), n=5))
## axis.bee.abund <-  axis.bee.abund
## HB abund
labs.HB.abund <- (pretty(c(0, spec.orig$Net_HBAbundance), n=5))
axis.HB.abund <-  standardize.axis(labs.HB.abund, spec.orig$Net_HBAbundance)
## bombus abund
labs.bombus.abund <- (pretty(c(0, spec.orig$Net_BombusAbundance), n=5))
axis.bombus.abund <-  standardize.axis(labs.bombus.abund, spec.orig$Net_BombusAbundance)
## non hb non bombus abund
labs.bee.abund <- (pretty(c(0, spec.orig$Net_NonBombusHBAbundance), n=5))
axis.bee.abund <-  standardize.axis(labs.bee.abund, spec.orig$Net_NonBombusHBAbundance)
## bee diversity
labs.bee.div <- (pretty(c(0, spec.orig$Net_BeeDiversity), n=5))
axis.bee.div <-  standardize.axis(labs.bee.div,
                                  spec.orig$Net_BeeDiversity)
## Load the model data
load(file="saved/communityFit.Rdata")
## We want the standarized data for the predictions (spec.   ) depending on the model
data.par <- spec.net[spec.net$WeightsPar == 1,]
## ***********************************************************************
## bee community diversity and latitude
## ***********************************************************************

## Community level visuals

newdata.lat <- crossing(Lat =
                          seq(min(data.par$Lat),
                              max(data.par$Lat),
                              length.out=10),
                        Net_BeeDiversity = mean(data.par$Net_BeeDiversity),
                        rare.degree = mean(data.par$rare.degree),
                        MeanITD = mean(data.par$MeanITD),
                        MeanFloralAbund= mean(data.par$MeanFloralAbundance),
                        MeanFloralDiversity= mean(data.par$MeanFloralDiversity),
                        Net_HBAbundance= mean(data.par$Net_HBAbundance),
                        Net_NonBombusHBAbundance = mean(data.par$Net_NonBombusHBAbundance),
                        Net_BombusAbundance = mean(data.par$Net_BombusAbundance),
                        Site = "JC", 
                        SRDoy = mean(data.par$SRDoy),
                        Year = "2017"
)

## predict values based on generated data and model parameters
pred_lat <- fit.community %>% 
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
  geom_point(data=data.par,
             aes(y= Net_BeeDiversity, x=Lat),
             color="grey40", cex=2)

ggsave(beediv_lat, file="figures/Lat_beediv.pdf",
       height=4, width=5)

################################################################################
## Plant community diversity and latitude
################################################################################

newdata.lat <- crossing(Lat =
                          seq(min(data.par$Lat),
                              max(data.par$Lat),
                              length.out=10),
                        Net_BeeDiversity = mean(data.par$Net_BeeDiversity),
                        rare.degree = mean(data.par$rare.degree),
                        MeanITD = mean(data.par$MeanITD),
                        MeanFloralAbund= mean(data.par$MeanFloralAbundance),
                        MeanFloralDiversity= mean(data.par$MeanFloralDiversity),
                        Net_HBAbundance= mean(data.par$Net_HBAbundance),
                        Net_NonBombusHBAbundance = mean(data.par$Net_NonBombusHBAbundance),
                        Net_BombusAbundance = mean(data.par$Net_BombusAbundance),
                        Site = "JC", 
                        SRDoy = mean(data.par$SRDoy),
                        Year = "2017"
)

## predict values based on generated data and model parameters
pred_lat <- fit.community %>% 
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
  geom_point(data=data.par,
             aes(y= MeanFloralDiversity, x=Lat),
             color="grey40", cex=2)

ggsave(plantdiv_lat, file="figures/Lat_floraldiv.pdf",
       height=4, width=5)
################################################################################
## Plant diversity and bee diversity
################################################################################
newdata.div <- crossing(Net_BeeDiversity =
                          seq(min(data.par$Net_BeeDiversity),
                              max(data.par$Net_BeeDiversity),
                              length.out=10),
                        Lat = mean(data.par$Lat),
                        rare.degree = mean(data.par$rare.degree),
                        MeanITD = mean(data.par$MeanITD),
                        MeanFloralAbund= mean(data.par$MeanFloralAbundance),
                        MeanFloralDiversity= mean(data.par$MeanFloralDiversity),
                        Net_HBAbundance= mean(data.par$Net_HBAbundance),
                        Net_NonBombusHBAbundance = mean(data.par$Net_NonBombusHBAbundance),
                        Net_BombusAbundance = mean(data.par$Net_BombusAbundance),
                        Site = "JC", 
                        SRDoy = mean(data.par$SRDoy),
                        Year = "2017"
)

## predict values based on generated data and model parameters
pred_div <- fit.community %>% 
  epred_draws(newdata = newdata.div,
              resp = "MeanFloralDiversity")

## to see range of predicted values
pred_div %>%
  group_by(MeanFloralDiversity) %>%
  summarise(mean(.epred))

plantdiv_beediv <- ggplot(pred_div, aes(x = .epred, y = MeanFloralDiversity)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Bee Species Diversity", y = "Flowering Species Diversity",
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
  geom_point(data=data.par,
             aes(y= MeanFloralDiversity, x=Lat),
             color="grey40", cex=2)

ggsave(plantdiv_beediv, file="figures/Beediv_floraldiv.pdf",
       height=4, width=5)

################################################################################
## Bee abundance and lat
################################################################################
newdata.lat <- crossing(Lat =
                          seq(min(data.par$Lat),
                              max(data.par$Lat),
                              length.out=10),
                        Net_BeeDiversity = mean(data.par$Net_BeeDiversity),
                        rare.degree = mean(data.par$rare.degree),
                        MeanITD = mean(data.par$MeanITD),
                        MeanFloralAbundance= mean(data.par$MeanFloralAbundance),
                        MeanFloralDiversity= mean(data.par$MeanFloralDiversity),
                        Net_HBAbundance= mean(data.par$Net_HBAbundance),
                        Net_NonBombusHBAbundance = mean(data.par$Net_NonBombusHBAbundance),
                        Net_BombusAbundance = mean(data.par$Net_BombusAbundance),
                        Site = "JC", 
                        SRDoy = mean(data.par$SRDoy),
                        Year = "2017"
)

## predict values based on generated data and model parameters
pred_lat <- fit.community %>% 
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
  geom_point(data=data.par,
             aes(y= Net_BombusAbundance, x=Lat),
             color="grey40", cex=2)

ggsave(bombusabun_lat, file="figures/Lat_bombus_abudance.pdf",
       height=4, width=5)

## Honeybee abundance
newdata.lat <- crossing(Lat =
                          seq(min(data.par$Lat),
                              max(data.par$Lat),
                              length.out=10),
                        Net_BeeDiversity = mean(data.par$Net_BeeDiversity),
                        rare.degree = mean(data.par$rare.degree),
                        MeanITD = mean(data.par$MeanITD),
                        MeanFloralAbundance= mean(data.par$MeanFloralAbundance),
                        MeanFloralDiversity= mean(data.par$MeanFloralDiversity),
                        Net_HBAbundance= mean(data.par$Net_HBAbundance),
                        Net_NonBombusHBAbundance = mean(data.par$Net_NonBombusHBAbundance),
                        Net_BombusAbundance = mean(data.par$Net_BombusAbundance),
                        Site = "JC", 
                        SRDoy = mean(data.par$SRDoy),
                        Year = "2017"
)

## predict values based on generated data and model parameters
pred_lat <- fit.community %>% 
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
  geom_point(data=data.par,
             aes(y= Net_HBAbundance, x=Lat),
             color="grey40", cex=2)

ggsave(HBabun_lat, file="figures/Lat_HB_abudance.pdf",
       height=4, width=5)