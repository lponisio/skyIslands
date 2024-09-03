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
bombus.par <- spec.bombus[spec.bombus$WeightsSp == 1,]
data.site <- spec.net[spec.net$Weights == 1,]

## bombus.par <- spec.bombus[spec.bombus$WeightsPar == 1,]

## https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/

## ***************************************************************************
## Crithidia ~ bee diversity
newdata.beediv <- crossing(Net_BeeDiversity =
                               seq(min(bombus.par$Net_BeeDiversity),
                                   max(bombus.par$Net_BeeDiversity),
                                   length.out=10),
                           rare.degree = 0,
                           Net_BeeAbundance = 0,
                           MeanFloralAbundance = 0,
                           MeanFloralDiversity = 0,
                           Lat = 0,
                           Site = "SC", 
                           GenusSpecies = "Bombus centralis",
                           WeightsPar=1
                           )

## predict values based on generated data and model parameters
pred_beediv <- fit.parasite %>% 
  epred_draws(newdata = newdata.beediv,
              resp = "CrithidiaPresence")

pred_beediv<- pred_beediv %>% mutate(bee = "Bombus")
## to see range of predicted values
pred_beediv %>%
  group_by(Net_BeeDiversity) %>%
  summarise(mean(.epred))

pred_beediv <- rbind(pred_beediv, pred_beediv2)

p1.parasite <- ggplot(pred_beediv, aes(x = Net_BeeDiversity, y = .epred, color = bee)) +
  stat_lineribbon(aes(linetype = bee), alpha = .4) +
  scale_fill_brewer(palette = "OrRd") +
  scale_color_manual(values = c("black", "black")) +
  labs(x = "Bee community diversity", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.bee.div,
    labels =  labs.bee.div) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16))

    ##theme_dark_black()+
    geom_jitter(data=bombus.par,
            aes(y= SpCrithidiaParasitismRate, x=Net_BeeDiversity,
                colour = GenusSpecies, #shape = Site,
                size=SpScreened), width=0.05) +
    guides(size = guide_legend(title = "Individuals screened"), 
           colour = guide_legend(title = ""),fill = "none")
    
ggsave(p1.parasite, file="figures/parasite_beediv.pdf",
           height=5, width=10)
## ***************************************************************************
## crithidia ~ bumble bee abundance
newdata.bombusabund <- crossing(Net_BeeAbundance =
                                    seq(min(bombus.par$Net_BeeAbundance),
                                        max(bombus.par$Net_BeeAbundance),
                                        length.out=10),
                                Lat = 0,
                                rare.degree = 0,
                                Net_BeeDiversity =0,
                                MeanFloralAbundance = 0,
                                MeanFloralDiversity = 0,
                                Site = "JC", 
                                GenusSpecies = "Bombus centralis",
                                WeightsPar=1
                                )

## predict values based on generated data and model parameters
pred_bombusabund <- fit.parasite %>% 
  epred_draws(newdata = newdata.bombusabund,
              resp = "CrithidiaPresence")

## to see range of predicted values
pred_bombusabund %>%
  group_by(Net_BeeAbundance) %>%
  summarise(mean(.epred))

p2.parasite <- ggplot(pred_bombusabund, aes(x = Net_BeeAbundance, y = .epred)) +
    stat_lineribbon() +
    scale_fill_brewer(palette = "Blues") +
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
    geom_jitter(data=bombus.par,
                aes(y= SpCrithidiaParasitismRate, x=Net_BeeAbundance,
                    colour = GenusSpecies, #shape = Site,
                    size=SpScreened), width=0.05) +
  guides(size = guide_legend(title = "Individuals screened"), 
         colour = guide_legend(title = ""),fill = "none")

parasite.comm <- ggarrange(p1.parasite, p2.parasite, nrow=1,
                          common.legend = TRUE,
                          legend="bottom")
ggsave(parasite.comm, file="figures/parasite_beecomm.pdf",
       height=5, width=10)


## ***************************************************************************
## Crithidia ~ floral diversity
newdata.floraldiv <- crossing(MeanFloralDiversity =
                               seq(min(bombus.par$MeanFloralDiversity),
                                   max(bombus.par$MeanFloralDiversity),
                                   length.out=10),
                           rare.degree = 0,
                           Lat = 0,
                           Net_BeeAbundance = 0,
                           MeanFloralAbundance = 0,
                           Net_BeeDiversity= 0,
                           Site = "SC", 
                           GenusSpecies = "Bombus centralis",
                           WeightsPar=1
                           )

## predict values based on generated data and model parameters
pred_floraldiv <- fit.parasite %>% 
  epred_draws(newdata = newdata.floraldiv,
              resp = "CrithidiaPresence")

## to see range of predicted values
pred_floraldiv %>%
  group_by(MeanFloralDiversity) %>%
  summarise(mean(.epred))

p3.parasite <- ggplot(pred_floraldiv, aes(x = MeanFloralDiversity, y = .epred)) +
    stat_lineribbon() +
    scale_fill_brewer(palette = "Blues") +
    labs(x = "Floral community diversity", y = "Crithidia prevalence",
         fill = "Credible interval") +
    theme_ms() +
    theme(legend.position = "bottom") +
    scale_x_continuous(
        breaks = axis.flower.div,
        labels =  labs.flower.div) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
    ##theme_dark_black()+
    geom_jitter(data=bombus.par,
            aes(y= SpCrithidiaParasitismRate, x=MeanFloralDiversity,
                colour = GenusSpecies, #shape = Site,
                size=SpScreened), width=0.05) +
  guides(size = guide_legend(title = "Individuals screened"), 
         colour = guide_legend(title = ""),fill = "none")

parasite.comm2 <- ggarrange(p1.parasite, p3.parasite, nrow=1,
                          common.legend = TRUE,
                          legend="bottom")
ggsave(parasite.comm2, file="figures/parasite_diversity.pdf",
       height=5, width=10)


## ***************************************************************************
## crithidia ~ diet breadth 
## ***************************************************************************

## ***************************************************************************
## crithidia ~ degree
newdata.degree <- crossing(rare.degree =
                             seq(min(bombus.par$rare.degree),
                                 max(bombus.par$rare.degree),
                                 length.out=10),
                           Lat = 0,
                           Net_BeeDiversity= 0,
                           Net_BeeAbundance = 0, 
                           MeanFloralAbundance = 0,
                           MeanFloralDiversity = 0,
                           Site = "JC", 
                           GenusSpecies = "Bombus centralis",
                           WeightsPar=1
)

## predict values based on generated data and model parameters
pred_degree <- fit.parasite %>% 
  epred_draws(newdata = newdata.degree,
              resp = "CrithidiaPresence")

## to see range of predicted values
pred_degree %>%
  group_by(rare.degree) %>%
  summarise(mean(.epred))

p4.parasite <- ggplot(pred_degree, aes(x = rare.degree, y = .epred)) +
    stat_lineribbon() +
    scale_fill_brewer(palette = "Blues") +
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
    geom_jitter(data=bombus.par,
                aes(y= SpCrithidiaParasitismRate, x=rare.degree,
                    colour = GenusSpecies, #shape = Site,
                    size=SpScreened), width=0.05) +
  guides(size = guide_legend(title = "Individuals screened"), 
         colour = guide_legend(title = ""),fill = "none")

parasite.traits <- ggarrange(p2.parasite, p4.parasite, nrow=1,
                          common.legend = TRUE,
                          legend="bottom")
ggsave(parasite.traits, file="figures/parasite_traits.pdf",
       height=5, width=10)
## ***********************************************************************
## lat and crithidia
## ***********************************************************************

newdata.lat <- crossing(Lat =
                          seq(min(bombus.par$Lat),
                              max(bombus.par$Lat),
                              length.out=10),
                        rare.degree = 0,
                        Net_BeeDiversity= 0,
                        Net_BeeAbundance = 0, 
                        MeanFloralAbundance = 0,
                        MeanFloralDiversity = 0,
                        Site = "JC", 
                        GenusSpecies = "Bombus centralis",
                        WeightsPar=1
)

## predict values based on generated data and model parameters
pred_lat <- fit.parasite %>% 
  epred_draws(newdata = newdata.lat,
              resp = "CrithidiaPresence")

## to see range of predicted values
pred_lat %>%
  group_by(Lat) %>%
  summarise(mean(.epred))

p5.parasite <- ggplot(pred_lat, aes(x = Lat, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
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
  geom_point(data=bombus.par,
              aes(y= SpCrithidiaParasitismRate, x= Lat),  
             color="grey40", cex=2)



ggsave(p5.parasite, file="figures/Lat_crithidia.pdf",
       height=4, width=5)

# ****************************************************************************
#Apicystis spp
## ***************************************************************************
## Apicystis ~ bee diversity
newdata.beediv <- crossing(Net_BeeDiversity =
                             seq(min(bombus.par$Net_BeeDiversity),
                                 max(bombus.par$Net_BeeDiversity),
                                 length.out=10),
                           rare.degree = 0,
                           Net_BeeAbundance = 0,
                           MeanFloralAbundance = 0,
                           MeanFloralDiversity = 0,
                           Lat = 0,
                           Site = "SC", 
                           GenusSpecies = "Bombus centralis",
                           WeightsPar=1
)

## predict values based on generated data and model parameters
pred_beediv <- fit.parasite %>% 
  epred_draws(newdata = newdata.beediv,
              resp = "ApicystisSpp")

## to see range of predicted values
pred_beediv %>%
  group_by(Net_BeeDiversity) %>%
  summarise(mean(.epred))

p6.parasite <- ggplot(pred_beediv, aes(x = Net_BeeDiversity, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Bee community diversity", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.bee.div,
    labels =  labs.bee.div) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  ##theme_dark_black()+
  geom_jitter(data=bombus.par,
              aes(y= SpApicystisParasitismRate, x=Net_BeeDiversity,
                  colour = GenusSpecies, #shape = Site,
                  size=SpScreened), width=0.05) +
  guides(size = guide_legend(title = "Individuals screened"), 
         colour = guide_legend(title = ""),fill = "none")

## ***************************************************************************
## Apicystis ~ bumble bee abundance
newdata.bombusabund <- crossing(Net_BeeAbundance =
                                  seq(min(bombus.par$Net_BeeAbundance),
                                      max(bombus.par$Net_BeeAbundance),
                                      length.out=10),
                                Lat = 0,
                                rare.degree = 0,
                                Net_BeeDiversity =0,
                                MeanFloralAbundance = 0,
                                MeanFloralDiversity = 0,
                                Site = "JC", 
                                GenusSpecies = "Bombus centralis",
                                WeightsPar=1
)

## predict values based on generated data and model parameters
pred_bombusabund <- fit.parasite %>% 
  epred_draws(newdata = newdata.bombusabund,
              resp = "ApicystisSpp")

## to see range of predicted values
pred_bombusabund %>%
  group_by(Net_BeeAbundance) %>%
  summarise(mean(.epred))

p7.parasite <- ggplot(pred_bombusabund, aes(x = Net_BeeAbundance, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Bee abundance (log)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.bee.abund,
    labels =  labs.bee.abund) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  geom_jitter(data=bombus.par,
              aes(y= SpApicystisParasitismRate, x=Net_BeeAbundance,
                  colour = GenusSpecies, #shape = Site,
                  size=SpScreened), width=0.05) +
  guides(size = guide_legend(title = "Individuals screened"), 
         colour = guide_legend(title = ""),fill = "none")

apicystis.comm <- ggarrange(p5.parasite, p6.parasite, nrow=1,
                           common.legend = TRUE,
                           legend="bottom")
ggsave(apicystis.comm, file="figures/Apicystis_beecomm.pdf",
       height=5, width=10)


## ***************************************************************************
## Apicystis ~ floral diversity
newdata.floraldiv <- crossing(MeanFloralDiversity =
                                seq(min(bombus.par$MeanFloralDiversity),
                                    max(bombus.par$MeanFloralDiversity),
                                    length.out=10),
                              rare.degree = 0,
                              Lat = 0,
                              Net_BeeAbundance = 0,
                              MeanFloralAbundance = 0,
                              Net_BeeDiversity= 0,
                              Site = "SC", 
                              GenusSpecies = "Bombus centralis",
                              WeightsPar=1
)

## predict values based on generated data and model parameters
pred_floraldiv <- fit.parasite %>% 
  epred_draws(newdata = newdata.floraldiv,
              resp = "ApicystisSpp")

## to see range of predicted values
pred_floraldiv %>%
  group_by(MeanFloralDiversity) %>%
  summarise(mean(.epred))

p8.parasite <- ggplot(pred_floraldiv, aes(x = MeanFloralDiversity, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Floral community diversity", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.flower.div,
    labels =  labs.flower.div) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  ##theme_dark_black()+
  geom_jitter(data=bombus.par,
              aes(y= SpApicystisParasitismRate, x=MeanFloralDiversity,
                  colour = GenusSpecies, #shape = Site,
                  size=SpScreened), width=0.05) +
  guides(size = guide_legend(title = "Individuals screened"), 
         colour = guide_legend(title = ""),fill = "none")

apicystis.comm2 <- ggarrange(p5.parasite, p7.parasite, nrow=1,
                            common.legend = TRUE,
                            legend="bottom")
ggsave(apicystis.comm2, file="figures/apicystis_diversity.pdf",
       height=5, width=10)


## ***************************************************************************
## Apicystis ~ diet breadth 
## ***************************************************************************

## ***************************************************************************
## Apicystis ~ degree
newdata.degree <- crossing(rare.degree =
                             seq(min(bombus.par$rare.degree),
                                 max(bombus.par$rare.degree),
                                 length.out=10),
                           Lat = 0,
                           Net_BeeDiversity= 0,
                           Net_BeeAbundance = 0, 
                           MeanFloralAbundance = 0,
                           MeanFloralDiversity = 0,
                           Site = "JC", 
                           GenusSpecies = "Bombus centralis",
                           WeightsPar=1
)

## predict values based on generated data and model parameters
pred_degree <- fit.parasite %>% 
  epred_draws(newdata = newdata.degree,
              resp = "ApicystisSpp")

## to see range of predicted values
pred_degree %>%
  group_by(rare.degree) %>%
  summarise(mean(.epred))

p9.parasite <- ggplot(pred_degree, aes(x = rare.degree, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
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
  geom_jitter(data=bombus.par,
              aes(y= SpApicystisParasitismRate, x=rare.degree,
                  colour = GenusSpecies, #shape = Site,
                  size=SpScreened), width=0.05) +
  guides(size = guide_legend(title = "Individuals screened"), 
         colour = guide_legend(title = ""),fill = "none")

apicystis.traits <- ggarrange(p6.parasite, p8.parasite, nrow=1,
                             common.legend = TRUE,
                             legend="bottom")

ggsave(apicystis.traits, file="figures/apicystis_traits.pdf",
       height=5, width=10)



#ggsave(p3.parasite, file="figures/parasite_Apicystisspp_beeDiv.jpg",
#       height=4, width=5)

#ggsave(p4.parasite, file="figures/parasite_Apicystisspp_bombusAbund.jpg",
#      height=4, width=5)


## ***********************************************************************
## lat and apicystis
## ***********************************************************************

newdata.lat <- crossing(Lat =
                             seq(min(bombus.par$Lat),
                                 max(bombus.par$Lat),
                                 length.out=10),
                           rare.degree = 0,
                           Net_BeeDiversity= 0,
                           Net_BeeAbundance = 0, 
                           MeanFloralAbundance = 0,
                           MeanFloralDiversity = 0,
                           Site = "JC", 
                           GenusSpecies = "Bombus centralis",
                           WeightsPar=1
)

## predict values based on generated data and model parameters
pred_lat <- fit.parasite %>% 
  epred_draws(newdata = newdata.lat,
              resp = "ApicystisSpp")

## to see range of predicted values
pred_lat %>%
  group_by(Lat) %>%
  summarise(mean(.epred))

p10.parasite <- ggplot(pred_lat, aes(x = Lat, y = .epred)) +
  stat_lineribbon() +
  scale_fill_manual(values = c("lightgoldenrod", "goldenrod1", "darkgoldenrod3")) +
  labs(x = "Latitude (log)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  geom_jitter(data=bombus.par,
              aes(y= SpApicystisParasitismRate, x= Lat), 
              color="grey40", cex=2) 



ggsave(p10.parasite, file="figures/Lat_apicystis.pdf",
       height=4, width=5)

lat_parasite<- ggarrange(p5.parasite, p11.parasite, p10.parasite, p12.parasite,#plots that are going to be included in this multipanel figure
                          labels = c("A", "C", "B","D"), #labels given each panel 
                          ncol = 2, nrow = 2)
 
ggsave(lat_parasite, file="figures/lat_parasites.pdf",
       height=8, width=12)


