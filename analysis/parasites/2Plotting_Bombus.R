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
## scaling/unscaling labs
## ***********************************************************************
spec.uni <- spec.orig[spec.orig$Weights ==1,]
## lat (logged)
labs.lat.x <- pretty(c(spec.uni$Lat),
                      n=10)
axis.lat.x <-  standardize.axis(labs.lat.x, spec.uni$Lat)

## bloom abundance (not logged)
labs.bloom.abund <- (pretty(c(spec.uni$MeanFloralAbundance), n=6))
axis.bloom.abund <-  standardize.axis(labs.bloom.abund,
                                      spec.uni$MeanFloralAbundance)
## flower div (not logged)
labs.flower.div <- (pretty(spec.uni$MeanFloralDiversity, n=5))
axis.flower.div <-  standardize.axis(labs.flower.div,
                                     spec.uni$MeanFloralDiversity)
## HB abund (logged + 1)
labs.HB.abund <- (pretty(c(spec.uni$Net_HBAbundance), n=5))
axis.HB.abund <-  standardize.axis(labs.HB.abund, spec.uni$Net_HBAbundance)
## bombus abund (logged + 1)
labs.bombus.abund <- (pretty(c(spec.uni$Net_BombusAbundance), n=5))
axis.bombus.abund <-  standardize.axis(labs.bombus.abund, spec.uni$Net_BombusAbundance)
## all bee abund (logged)
labs.bee.abund <- (pretty(c(spec.uni$Net_BeeAbundance), n=5))
axis.bee.abund <-  standardize.axis(labs.bee.abund, spec.uni$Net_BeeAbundance)
## bee diversity (not logged)
labs.bee.div <- (pretty(c(spec.uni$Net_BeeDiversity), n=5))
axis.bee.div <-  standardize.axis(labs.bee.div,
                                  spec.uni$Net_BeeDiversity)

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

load(file="../../../skyIslands_saved/parasite-results/saved/parasiteFit_apis_CrithidiaPresenceApicystisSpp_lat_ss.Rdata")
fit.apis <- fit.parasite
## https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/

## ***************************************************************************
## Crithidia ~ bee diversity
newdata.beediv <- crossing(Net_BeeDiversity =
                               seq(min(spec.uni$Net_BeeDiversity),
                                   max(spec.uni$Net_BeeDiversity),
                                   length.out=10),
                           rare.degree = 0,
                           Net_BeeAbundance = 0,
                           MeanFloralAbundance = 0,
                           MeanFloralDiversity = 0,
                           Lat = 0,
                           Site = "JC", 
                           GenusSpecies = "Bombus centralis",
                           WeightsPar=1
                           )

## predict values based on generated data and model parameters
pred_beediv <- fit.bombus %>% 
  epred_draws(newdata = newdata.beediv,
              resp = "CrithidiaPresence")

pred_beediv<- pred_beediv %>% mutate(bee = "Bombus")
## to see range of predicted values
pred_beediv %>%
  group_by(Net_BeeDiversity) %>%
  summarise(mean(.epred))

## Create the same fake data but using the Apis model

newdata.beediv2 <- crossing(Net_BeeDiversity =
                               seq(min(spec.uni$Net_BeeDiversity),
                                   max(spec.uni$Net_BeeDiversity),
                                   length.out=10),
                           rare.degree = 0,
                           Net_BombusAbundance = 0,
                           Net_HBAbundance = 0,
                           MeanFloralAbundance = 0,
                           MeanFloralDiversity = 0,
                           Lat = 0,
                           Site = "JC", 
                           GenusSpecies = "Apis mellifera",
                           WeightsPar=1
                           )

## predict values based on generated data and model parameters
pred_beediv2 <- fit.apis %>% 
  epred_draws(newdata = newdata.beediv2,
              resp = "CrithidiaPresence")
pred_beediv2<- pred_beediv2 %>% mutate(bee = "Apis")
## to see range of predicted values
pred_beediv2 %>%
  group_by(Net_BeeDiversity) %>%
  summarise(mean(.epred))

## Merge the two fake data
pred_beediv <- rbind(pred_beediv, pred_beediv2)

p1.parasite <- ggplot(pred_beediv, aes(x = Net_BeeDiversity, y = .epred, fill = bee)) +
  stat_lineribbon(aes(linetype = bee, color = bee), alpha = .6, .width = .95) +
  scale_fill_manual(values = c("darkgoldenrod3", "#3182bd"), labels = c("Apis 0.95", "Bombus 0.95")) +
  scale_color_manual(values = c("black", "black")) +
  labs(x = "Bee diversity", y = "Crithidia prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.bee.div,
    labels =  labs.bee.div) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        text = element_text(size=16))

    ##theme_dark_black()+
    # geom_jitter(data=spec.uni,
    #         aes(y= SpCrithidiaParasitismRate, x=Net_BeeDiversity,
    #             colour = GenusSpecies, #shape = Site,
    #             size=SpScreened), width=0.05) +
    # guides(size = guide_legend(title = "Individuals screened"), 
    #        colour = guide_legend(title = ""),fill = "none")
    
ggsave(p1.parasite, file="figures/parasite_beediv_Crithidia.pdf",
           height=5, width=10)

################################################################################
## Apicystis ~ bee diversity Bombus
################################################################################

newdata.beediv <- crossing(Net_BeeDiversity =
                             seq(min(spec.uni$Net_BeeDiversity),
                                 max(spec.uni$Net_BeeDiversity),
                                 length.out=10),
                           rare.degree = 0,
                           Net_BeeAbundance = 0,
                           MeanFloralAbundance = 0,
                           MeanFloralDiversity = 0,
                           Lat = 0,
                           Site = "JC", 
                           GenusSpecies = "Bombus centralis",
                           WeightsPar=1
)

## predict values based on generated data and model parameters
pred_beediv3 <- fit.bombus %>% 
  epred_draws(newdata = newdata.beediv,
              resp = "ApicystisSpp")
#pred_beediv3 <- pred_beediv3 %>% mutate(bee = "Bombus")
## to see range of predicted values
pred_beediv3 %>%
  group_by(Net_BeeDiversity) %>%
  summarise(mean(.epred))
## Load Apis and create the predicted data

# 
# newdata.beediv <- crossing(Net_BeeDiversity =
#                              seq(min(apis.par$Net_BeeDiversity),
#                                  max(apis.par$Net_BeeDiversity),
#                                  length.out=10),
#                            rare.degree = 0,
#                            Net_BombusAbundance = 0,
#                            Net_HBAbundance = 0,
#                            MeanFloralAbundance = 0,
#                            MeanFloralDiversity = 0,
#                            Lat = 0,
#                            Site = "JC", 
#                            GenusSpecies = "Bombus centralis",
#                            WeightsPar=1
# )
# 
# ## predict values based on generated data and model parameters
# pred_beediv4 <- fit.apis %>% 
#   epred_draws(newdata = newdata.beediv,
#               resp = "ApicystisSpp")
# pred_beediv4 <- pred_beediv4 %>% mutate(bee = "Apis")
# ## to see range of predicted values
# pred_beediv4 %>%
#   group_by(Net_BeeDiversity) %>%
#   summarise(mean(.epred))
# 
# ## Merge the two data 
# pred_beediv3 <- rbind(pred_beediv3, pred_beediv4)

p2.parasite <- ggplot(pred_beediv3, aes(x = Net_BeeDiversity, y = .epred)) +
  stat_lineribbon(.width = .95, linetype = "dashed", alpha = .6) +
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
        text = element_text(size=16))

ggsave(p2.parasite, file="figures/parasite_beediv_Apicystis.pdf",
       height=5, width=10)

## ***************************************************************************
## Crithidia ~ floral diversity

newdata.floraldiv <- crossing(MeanFloralDiversity =
                               seq(min(spec.uni$MeanFloralDiversity),
                                   max(spec.uni$MeanFloralDiversity),
                                   length.out=10),
                           rare.degree = 0,
                           Lat = 0,
                           Net_BeeAbundance = 0,
                           MeanFloralAbundance = 0,
                           Net_BeeDiversity= 0,
                           Site = "JC", 
                           GenusSpecies = "Bombus centralis",
                           WeightsPar=1
                           )

## predict values based on generated data and model parameters
pred_floraldiv <- fit.bombus %>% 
  epred_draws(newdata = newdata.floraldiv,
              resp = "CrithidiaPresence")
#pred_floraldiv<- pred_floraldiv %>% mutate(bee = "Bombus")
## to see range of predicted values
pred_floraldiv %>%
  group_by(MeanFloralDiversity) %>%
  summarise(mean(.epred))
## Create the same fake data but using the Apis model
## You will need to load the apis model

# 
# newdata.floraldiv <- crossing(MeanFloralDiversity =
#                                 seq(min(apis.par$MeanFloralDiversity),
#                                     max(apis.par$MeanFloralDiversity),
#                                     length.out=10),
#                               rare.degree = 0,
#                               Lat = 0,
#                               Net_BombusAbundance = 0,
#                               Net_HBAbundance = 0,
#                               MeanFloralAbundance = 0,
#                               Net_BeeDiversity= 0,
#                               Site = "JC", 
#                               GenusSpecies = "Bombus centralis",
#                               WeightsPar=1
# )
# 
# ## predict values based on generated data and model parameters
# pred_floraldiv2 <- fit.apis %>% 
#   epred_draws(newdata = newdata.floraldiv,
#               resp = "CrithidiaPresence")
# pred_floraldiv2 <- pred_floraldiv2 %>% mutate(bee = "Apis")
# ## to see range of predicted values
# pred_floraldiv2 %>%
#   group_by(MeanFloralDiversity) %>%
#   summarise(mean(.epred))
# 
# ## Merge the fake data of Bombus and Apis and use it to graph
# pred_floraldiv <- rbind(pred_floraldiv, pred_floraldiv2)

p3.parasite <- ggplot(pred_floraldiv, aes(x = MeanFloralDiversity, y = .epred)) +
  stat_lineribbon(.width = .95, linetype = "dashed", alpha = .6) +
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
          text = element_text(size=16)) 
    ##theme_dark_black()+

ggsave(p3.parasite, file="figures/parasite_floraldiv_Crithidia.pdf",
       height=5, width=10)

## ***************************************************************************
## Apicystis ~ floral diversity

newdata.floraldiv <- crossing(MeanFloralDiversity =
                                seq(min(spec.uni$MeanFloralDiversity),
                                    max(spec.uni$MeanFloralDiversity),
                                    length.out=10),
                              rare.degree = 0,
                              Lat = 0,
                              Net_BeeAbundance = 0,
                              MeanFloralAbundance = 0,
                              Net_BeeDiversity= 0,
                              Site = "JC", 
                              GenusSpecies = "Bombus centralis",
                              WeightsPar=1
)

## predict values based on generated data and model parameters
pred_floraldiv3 <- fit.bombus %>% 
  epred_draws(newdata = newdata.floraldiv,
              resp = "ApicystisSpp")
#pred_floraldiv3<- pred_floraldiv3 %>% mutate(bee = "Bombus")
## to see range of predicted values
pred_floraldiv3 %>%
  group_by(MeanFloralDiversity) %>%
  summarise(mean(.epred))

## Now do the same but with the apis model

# 
# newdata.floraldiv <- crossing(MeanFloralDiversity =
#                                 seq(min(apis.par$MeanFloralDiversity),
#                                     max(apis.par$MeanFloralDiversity),
#                                     length.out=10),
#                               rare.degree = 0,
#                               Lat = 0,
#                               Net_BombusAbundance = 0,
#                               Net_HBAbundance = 0,
#                               MeanFloralAbundance = 0,
#                               Net_BeeDiversity= 0,
#                               Site = "JC", 
#                               GenusSpecies = "Bombus centralis",
#                               WeightsPar=1
# )
# 
# ## predict values based on generated data and model parameters
# pred_floraldiv4 <- fit.apis %>% 
#   epred_draws(newdata = newdata.floraldiv,
#               resp = "ApicystisSpp")
# pred_floraldiv4<- pred_floraldiv4 %>% mutate(bee = "Apis")
# ## to see range of predicted values
# pred_floraldiv4 %>%
#   group_by(MeanFloralDiversity) %>%
#   summarise(mean(.epred))
# ## Merge the data 
# pred_floraldiv3 <- rbind(pred_floraldiv3, pred_floraldiv4)

p4.parasite <- ggplot(pred_floraldiv3, aes(x = MeanFloralDiversity, y = .epred)) +
  stat_lineribbon(.width = .95, linetype = "dashed", alpha = .6) +
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
        text = element_text(size=16))

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

newdata.bombusabund <- crossing(Net_BeeAbundance =
                                  seq(min(spec.uni$Net_BeeAbundance),
                                      max(spec.uni$Net_BeeAbundance),
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
pred_bombusabund <- fit.bombus %>% 
  epred_draws(newdata = newdata.bombusabund,
              resp = "CrithidiaPresence")

## to see range of predicted values
pred_bombusabund %>%
  group_by(Net_BeeAbundance) %>%
  summarise(mean(.epred))

p5.parasite <- ggplot(pred_bombusabund, aes(x = Net_BeeAbundance, y = .epred)) +
  stat_lineribbon(alpha = .6, .width = .95) +
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
        text = element_text(size=16)) 
  # geom_jitter(data=spec.uni,
  #             aes(y= SpCrithidiaParasitismRate, x=Net_BeeAbundance,
  #                 colour = GenusSpecies, #shape = Site,
  #                 size=SpScreened), width=0.05) +
  # guides(size = guide_legend(title = "Individuals screened"), 
  #        colour = guide_legend(title = ""),fill = "none")


ggsave(p5.parasite, file="figures/crithidia_beeabundance_bombus.pdf",
       height=4, width=6)
################################################################################
## Apicysits ~ bee abundance
################################################################################
newdata.bombusabund <- crossing(Net_BeeAbundance =
                                  seq(min(spec.uni$Net_BeeAbundance),
                                      max(spec.uni$Net_BeeAbundance),
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
pred_bombusabund <- fit.bombus %>% 
  epred_draws(newdata = newdata.bombusabund,
              resp = "ApicystisSpp")

## to see range of predicted values
pred_bombusabund %>%
  group_by(Net_BeeAbundance) %>%
  summarise(mean(.epred))

p6.parasite <- ggplot(pred_bombusabund, aes(x = Net_BeeAbundance, y = .epred)) +
  stat_lineribbon(.width = .95, alpha = 0.6) +
  scale_fill_manual(values = "#3182bd", labels ="Bombus 0.95") +
  labs(x = "Bee abundance (log)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.bee.abund,
    labels =  labs.bee.abund) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) 
  # geom_jitter(data=spec.uni,
  #             aes(y= SpApicystisParasitismRate, x=Net_BeeAbundance,
  #                 colour = GenusSpecies, #shape = Site,
  #                 size=SpScreened), width=0.05) +
  # guides(size = guide_legend(title = "Individuals screened"), 
  #        colour = guide_legend(title = ""),fill = "none")

ggsave(p6.parasite, file="figures/Apicystis_beeabun_bombus.pdf",
       height=4, width=6)

parasite.amplification <- ggarrange(p5.parasite, p6.parasite, nrow=1,
                                    labels = c("A", "B"),
                                    common.legend = TRUE,
                                    legend="bottom")

ggsave(parasite.amplification, file="figures/parasite_amplification.pdf",
       height=6, width=10)

## ***************************************************************************
## parasites ~ diet breadth 
## ***************************************************************************

## ***************************************************************************
## crithidia ~ degree
newdata.degree <- crossing(rare.degree =
                             seq(min(bombus.par$rare.degree, na.rm=TRUE),
                                 max(bombus.par$rare.degree, na.rm=TRUE),
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
pred_degree <- fit.bombus %>% 
  epred_draws(newdata = newdata.degree,
              resp = "CrithidiaPresence")

## to see range of predicted values
pred_degree %>%
  group_by(rare.degree) %>%
  summarise(mean(.epred))

p7.parasite <- ggplot(pred_degree, aes(x = rare.degree, y = .epred)) +
    stat_lineribbon(.width = 0.95, alpha = 0.6) +
  scale_fill_manual(values = "#3182bd", labels ="Bombus 0.95") +
    labs(x = "Degree", y = "Crithidia prevalence",
         fill = "Credible interval") +
    theme_ms() +
    theme(legend.position = "bottom") +
    scale_x_continuous(
        breaks = axis.degree,
        labels =  labs.degree) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) 
  #   geom_jitter(data=bombus.par,
  #               aes(y= SpCrithidiaParasitismRate, x=rare.degree,
  #                   colour = GenusSpecies, #shape = Site,
  #                   size=SpScreened), width=0.05) +
  # guides(size = guide_legend(title = "Individuals screened"), 
  #        colour = guide_legend(title = ""),fill = "none")


ggsave(p7.parasite, file="figures/crithidia_degree.pdf",
       height=5, width=10)

## Apicystis ~ degree
newdata.degree <- crossing(rare.degree =
                             seq(min(bombus.par$rare.degree, na.rm=TRUE),
                                 max(bombus.par$rare.degree, na.rm=TRUE),
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
pred_degree <- fit.bombus %>% 
  epred_draws(newdata = newdata.degree,
              resp = "ApicystisSpp")

## to see range of predicted values
pred_degree %>%
  group_by(rare.degree) %>%
  summarise(mean(.epred))

p8.parasite <- ggplot(pred_degree, aes(x = rare.degree, y = .epred)) +
  stat_lineribbon(.width = 0.95, alpha = 0.6) +
  scale_fill_manual(values = "#3182bd", labels ="Bombus 0.95") +
  labs(x = "Degree", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.degree,
    labels =  labs.degree) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) 
  # geom_jitter(data=bombus.par,
  #             aes(y= SpApicystisParasitismRate, x=rare.degree,
  #                 colour = GenusSpecies, #shape = Site,
  #                 size=SpScreened), width=0.05) +
  # guides(size = guide_legend(title = "Individuals screened"), 
  #        colour = guide_legend(title = ""),fill = "none")

ggsave(p8.parasite, file="figures/apicystis_degree.pdf",
       height=5, width=10)

parasite.traits <- ggarrange(p7.parasite, p8.parasite, nrow=1,
                             labels = c("A", "B"),
                             common.legend = TRUE,
                             legend="bottom")
ggsave(parasite.traits, file="figures/parasite_traits.pdf",
       height= 6, width=10)
## ***********************************************************************
## lat and crithidia
## ***********************************************************************

newdata.lat <- crossing(Lat =
                          seq(min(spec.uni$Lat),
                              max(spec.uni$Lat),
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
pred_lat <- fit.bombus %>% 
  epred_draws(newdata = newdata.lat,
              resp = "CrithidiaPresence")

## to see range of predicted values
pred_lat %>%
  group_by(Lat) %>%
  summarise(mean(.epred))

p9.parasite <- ggplot(pred_lat, aes(x = Lat, y = .epred)) +
  stat_lineribbon(.width = 0.95, alpha = 0.6) +
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
        plot.title = element_text(color = "black")) 
  # geom_point(data=bombus.par,
  #             aes(y= SpCrithidiaParasitismRate, x= Lat),  
  #            color="grey40", cex=2)



ggsave(p9.parasite, file="figures/Lat_crithidia.pdf",
       height=4, width=5)

## ***********************************************************************
## lat and apicystis
## ***********************************************************************

newdata.lat <- crossing(Lat =
                             seq(min(spec.uni$Lat),
                                 max(spec.uni$Lat),
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
pred_lat <- fit.bombus %>% 
  epred_draws(newdata = newdata.lat,
              resp = "ApicystisSpp")

## to see range of predicted values
pred_lat %>%
  group_by(Lat) %>%
  summarise(mean(.epred))

p10.parasite <- ggplot(pred_lat, aes(x = Lat, y = .epred)) +
  stat_lineribbon(.width = 0.95, alpha = 0.6) +
  scale_fill_manual(values = "darkgoldenrod3", labels ="Bombus 0.95") +
  labs(x = "Latitude (log)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) 
  # geom_jitter(data=bombus.par,
  #             aes(y= SpApicystisParasitismRate, x= Lat), 
  #             color="grey40", cex=2) 



ggsave(p10.parasite, file="figures/Lat_apicystis.pdf",
       height=4, width=5)

################################################################################
## crithida ~ lat in Apis
newdata.lat <- crossing(Lat =
                          seq(min(spec.uni$Lat),
                              max(spec.uni$Lat),
                              length.out=10),
                        rare.degree = 0,
                        Net_BeeDiversity= 0, 
                        Net_BombusAbundance = 0,
                        Net_HBAbundance = 0,
                        MeanFloralAbundance = 0,
                        MeanFloralDiversity = 0,
                        Site = "JC", 
                        GenusSpecies = "Bombus centralis",
                        WeightsPar=1
)

## predict values based on generated data and model parameters
pred_lat <- fit.apis %>% 
  epred_draws(newdata = newdata.lat,
              resp = "CrithidiaPresence")

## to see range of predicted values
pred_lat %>%
  group_by(Lat) %>%
  summarise(mean(.epred))

p11.parasite <- ggplot(pred_lat, aes(x = Lat, y = .epred)) +
  stat_lineribbon(.width = 0.95, alpha = 0.6) +
  scale_fill_manual(values = "#e6550d", labels ="Apis 0.95") +
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
        plot.title = element_text(color = "black")) 
# geom_jitter(data=bombus.par,
#             aes(y= SpCrithidiaParasitismRate, x= Lat), 
#             color="grey40", cex=2) 



ggsave(p11.parasite, file="figures/Lat_crithidia_apis.pdf",
       height=4, width=5)

################################################################################
## Lat and apicystis in apis
################################################################################
newdata.lat <- crossing(Lat =
                          seq(min(spec.uni$Lat),
                              max(spec.uni$Lat),
                              length.out=10),
                        rare.degree = 0,
                        Net_BeeDiversity= 0,
                        Net_BombusAbundance = 0,
                        Net_HBAbundance = 0,
                        MeanFloralAbundance = 0,
                        MeanFloralDiversity = 0,
                        Site = "JC", 
                        GenusSpecies = "Bombus centralis",
                        WeightsPar=1
)

## predict values based on generated data and model parameters
pred_lat <- fit.apis %>% 
  epred_draws(newdata = newdata.lat,
              resp = "ApicystisSpp")

## to see range of predicted values
pred_lat %>%
  group_by(Lat) %>%
  summarise(mean(.epred))

p12.parasite <- ggplot(pred_lat, aes(x = Lat, y = .epred)) +
  stat_lineribbon(.width = 0.95, alpha = 0.6) +
  scale_fill_manual(values = "#e6550d", labels ="Apis 0.95") +
  labs(x = "Latitude (log)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) 
# geom_jitter(data=bombus.par,
#             aes(y= SpApicystisParasitismRate, x= Lat), 
#             color="grey40", cex=2) 



ggsave(p12.parasite, file="figures/Lat_apicystis_apis.pdf",
       height=4, width=5)



lat_parasite<- ggarrange(p9.parasite, p11.parasite, p10.parasite, p12.parasite,#plots that are going to be included in this multipanel figure
                          labels = c("A", "C", "B","D"), #labels given each panel 
                          ncol = 2, nrow = 2)
 
ggsave(lat_parasite, file="figures/lat_parasites.pdf",
       height=8, width=12)


