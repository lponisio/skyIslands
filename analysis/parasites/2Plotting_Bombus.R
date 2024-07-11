rm(list=ls())
## setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')
## setwd("C:/Users/na_ma/Dropbox (University of Oregon)/skyIslands")
setwd('~/Dropbox (University of Oregon)/skyislands')
## Script for plotting all of the important explanatory variables.
library(ggpubr)

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
load(file="saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat.Rdata")

# We want the standarized data for the predictions (spec.data)
spec.bombus <- spec.net[spec.net$Genus == "Bombus",]
data.par <- spec.bombus[spec.bombus$WeightsSp == 1,]
data.site <- spec.net[spec.net$Weights == 1,]

## data.par <- spec.bombus[spec.bombus$WeightsPar == 1,]

## https://www.rensvandeschoot.com/tutorials/generalised-linear-models-with-brms/

## ***************************************************************************
## Crithidia ~ bee diversity
newdata.beediv <- crossing(Net_BeeDiversity =
                               seq(min(data.par$Net_BeeDiversity),
                                   max(data.par$Net_BeeDiversity),
                                   length.out=10),
                           rare.degree = 0,
                           MeanITD =0,
                           Net_BombusAbundance = 0,
                           MeanFloralAbundance = 0,
                           MeanFloralDiversity = 0,
                           Site = "JC", 
                           GenusSpecies = "Bombus centralis",
                           WeightsPar=1
                           )

## predict values based on generated data and model parameters
pred_beediv <- fit.parasite %>% 
  epred_draws(newdata = newdata.beediv,
              resp = "CrithidiaPresence")

## to see range of predicted values
pred_beediv %>%
  group_by(Net_BeeDiversity) %>%
  summarise(mean(.epred))

p1.parasite <- ggplot(pred_beediv, aes(x = Net_BeeDiversity, y = .epred)) +
    stat_lineribbon() +
    scale_fill_brewer(palette = "Blues") +
    labs(x = "Bee community diversity", y = "Crithidia Prevalence",
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
    geom_jitter(data=data.par,
                aes(y= SpCrithidiaParasitismRate, x=Net_BeeDiversity,
                    colour = GenusSpecies, #shape = Site,
                    size=SpScreened), width=0.05) +
    guides(
        size = guide_legend(
            title = "Individuals screened",
            ),
        colour = guide_legend(
            title = "Species"
        ),
        fill = "none"
    )
#+ facet_wrap(~Year)

## ***************************************************************************
## crithidia ~ bumble bee abundance
newdata.bombusabund <- crossing(Net_BombusAbundance =
                                    seq(min(data.par$Net_BombusAbundance),
                                        max(data.par$Net_BombusAbundance),
                                        length.out=10),
                                Lat = 0,
                                rare.degree = 0,
                                MeanITD = 0,
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
  group_by(Net_BombusAbundance) %>%
  summarise(mean(.epred))

p2.parasite <- ggplot(pred_bombusabund, aes(x = Net_BombusAbundance, y = .epred)) +
    stat_lineribbon() +
    scale_fill_brewer(palette = "Blues") +
    labs(x = "Bombus abundance (log)", y = "Parasite Prevalence",
         fill = "Credible interval") +
    theme_ms() +
    theme(legend.position = "bottom") +
    scale_x_continuous(
        breaks = axis.bombus.abund,
        labels =  labs.bombus.abund) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
    geom_jitter(data=data.par,
                aes(y= SpCrithidiaParasitismRate, x=Net_BombusAbundance,
                    colour = GenusSpecies, #shape = Site,
                    size=SpScreened), width=0.05) +
    guides(
        size = guide_legend(
            title = "Individuals screened",
            ),
        colour = guide_legend(
            title = "Species"
        ),
        fill = "none"
    )

parasite.comm <- ggarrange(p1.parasite, p2.parasite, nrow=1,
                          common.legend = TRUE,
                          legend="bottom")
ggsave(parasite.comm, file="figures/parasite_beecomm.pdf",
       height=5, width=10)

## ***************************************************************************
## crithidia ~ bee traits 
## ***************************************************************************
## parasitism ~ meanitd
newdata.itd <- crossing(MeanITD =
                             seq(min(data.par$MeanITD),
                                 max(data.par$MeanITD),
                                 length.out=10),
                           rare.degree = 0,
                           Net_BeeDiversity= 0,
                           Net_BombusAbundance = 0, 
                           MeanFloralAbundance = 0,
                           MeanFloralDiversity = 0,
                           Site = "JC", 
                           GenusSpecies = "Bombus centralis",
                           WeightsPar=1
)

## predict values based on generated data and model parameters
pred_itd <- fit.parasite %>% 
  epred_draws(newdata = newdata.itd,
              resp = "CrithidiaPresence")

## to see range of predicted values
pred_itd %>%
  group_by(MeanITD) %>%
  summarise(mean(.epred))

p3.parasite <- ggplot(pred_itd, aes(x = MeanITD, y = .epred)) +
    stat_lineribbon() +
    scale_fill_brewer(palette = "Blues") +
    labs(x = "Body Size (ITD mm)", y = "Crithidia prevalence",
         fill = "Credible interval") +
    theme_ms() +
    theme(legend.position = "bottom") +
    scale_x_continuous(
        breaks = axis.itd,
        labels =  labs.itd) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
    ##theme_dark_black()+
    geom_jitter(data=data.par,
                aes(y= SpCrithidiaParasitismRate, x=MeanITD,
                    colour = GenusSpecies, #shape = Site,
                    size=SpScreened), width=0.05) +
    guides(
        size = guide_legend(
            title = "Individuals screened",
            ),
        colour = guide_legend(
            title = "Species"
        ),
        fill = "none"
    )
#+ facet_wrap(~Year)

## ***************************************************************************
## crithidia ~ degree
newdata.degree <- crossing(rare.degree =
                             seq(min(data.par$rare.degree),
                                 max(data.par$rare.degree),
                                 length.out=10),
                           MeanITD = 0,
                           Net_BeeDiversity= 0,
                           Net_BombusAbundance = 0, 
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
        breaks = axis.bombus.abund,
        labels =  labs.bombus.abund) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
    geom_jitter(data=data.par,
                aes(y= SpCrithidiaParasitismRate, x=rare.degree,
                    colour = GenusSpecies, #shape = Site,
                    size=SpScreened), width=0.05) +
    guides(
        size = guide_legend(
            title = "Individuals screened",
            ),
        colour = guide_legend(
            title = "Species"
        ),
        fill = "none"
    )

parasite.traits <- ggarrange(p3.parasite, p4.parasite, nrow=1,
                          common.legend = TRUE,
                          legend="bottom")
ggsave(parasite.traits, file="figures/parasite_traits.pdf",
       height=5, width=10)

# *******************************************************************

#Apicystis spp
newdata.beediv <- crossing(Net_BeeDiversity =
                             seq(min(data.par$Net_BeeDiversity),
                                 max(data.par$Net_BeeDiversity),
                                 length.out=10),
                           Lat = mean(data.par$Lat),
                           rare.degree = mean(data.par$rare.degree),
                           MeanITD = mean(data.par$MeanITD),
                           Net_BombusAbundance = mean(data.par$Net_BombusAbundance),
                           Site = "JC", 
                           GenusSpecies = "Bombus centralis"
)

## predict values based on generated data and model parameters
pred_beediv <- fit.parasite %>% 
  epred_draws(newdata = newdata.beediv,
              resp = "ApicystisSpp")

## to see range of predicted values
pred_beediv %>%
  group_by(Net_BeeDiversity) %>%
  summarise(mean(.epred))

p3.parasite <- ggplot(pred_beediv, aes(x = Net_BeeDiversity, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Bee community diversity", y = "Apicystis Spp Prevalence",
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
  geom_point(data=data.par,
             aes(y= SpApicystisParasitismRate, x=Net_BeeDiversity, colour = GenusSpecies),
             cex=2)

## parasitism ~ bumble bee abundance

newdata.bombusabund <- crossing(Net_BombusAbundance =
                                  seq(min(data.par$Net_BombusAbundance),
                                      max(data.par$Net_BombusAbundance),
                                      length.out=10),
                                Lat = mean(data.par$Lat),
                                rare.degree = mean(data.par$rare.degree),
                                MeanITD = mean(data.par$MeanITD),
                                Net_BeeDiversity = mean(data.par$Net_BeeDiversity),
                                Site = "JC", 
                                GenusSpecies = "Bombus centralis"
)

## predict values based on generated data and model parameters
pred_bombusabund <- fit.parasite %>% 
  epred_draws(newdata = newdata.bombusabund,
              resp = "ApicystisSpp")

## to see range of predicted values
pred_bombusabund %>%
  group_by(Net_BombusAbundance) %>%
  summarise(mean(.epred))

p4.parasite <- ggplot(pred_bombusabund, aes(x = Net_BombusAbundance, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Bumble bee Abundance", y = "Apicystis Spp Prevalence",
       fill = "Credible interval") +
  theme_ms() +
  theme(legend.position = "bottom") +
  #scale_x_continuous(
  #breaks = axis.bombus.abund,
  #labels =  labs.bombus.abund) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  ##theme_dark_black()+
  geom_point(data=data.par,
             aes(y= SpApicystisParasitismRate , x=Net_BombusAbundance),
             color="grey40", cex=2)



ggsave(p3.parasite, file="figures/parasite_Apicystisspp_beeDiv.jpg",
       height=4, width=5)

ggsave(p4.parasite, file="figures/parasite_Apicystisspp_bombusAbund.jpg",
       height=4, width=5)


## ***********************************************************************
## lat and bee diversity
## ***********************************************************************

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
                           GenusSpecies = "Bombus huntii",
                        SRDoy = mean(data.par$SRDoy),
                        Year = mean(data.par$Year)
)

## predict values based on generated data and model parameters
pred_lat <- fit.parasite %>% 
  epred_draws(newdata = newdata.lat,
              resp = "NetBeeDiversity")

## to see range of predicted values
pred_lat %>%
  group_by(Net_BeeDiversity) %>%
  summarise(mean(.epred))

p3.parasite <- ggplot(pred_lat, aes(x = .epred, y = Net_BeeDiversity)) +
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

ggsave(p3.parasite, file="figures/Lat_beediv.pdf",
       height=4, width=5)

## ***********************************************************************
## traits and parasites
## ***********************************************************************




ggsave(p3.parasite, file="figures/parasite_meanitd.pdf",
       height=4, width=5)

ggsave(p4.parasite, file="figures/parasite_rdegree.pdf",
       height=4, width=5)

parasite.all <- grid.arrange(p3.parasite, p4.parasite, ncol=2)

ggsave(parasite.all, file="figures/traits_parasite.pdf",
       height=4, width=10)


## ***********************************************************************
## flowers and parasites
## ***********************************************************************


p5.parasite  <- fit %>%
    spread_draws(b_MeanFloralDiversity_Intercept,
                 b_ParasitePresence_Intercept,
                 b_ParasitePresence_MeanFloralDiversity) %>%
    mutate(MeanFloralDiversity =
               list(seq(min(spec.all$MeanFloralDiversity[
                                           spec.all$WeightsPar == 1]),
                        max(spec.all$MeanFloralDiversity[
                                           spec.all$WeightsPar == 1]),
                        0.1))) %>%
    unnest(MeanFloralDiversity) %>%
    mutate(pred = exp(b_MeanFloralDiversity_Intercept +
                      b_ParasitePresence_Intercept +
                      b_ParasitePresence_MeanFloralDiversity*MeanFloralDiversity)/
               (1+exp(b_MeanFloralDiversity_Intercept +
                      b_ParasitePresence_Intercept +
                      b_ParasitePresence_MeanFloralDiversity*MeanFloralDiversity))) %>%
    group_by(MeanFloralDiversity) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = MeanFloralDiversity, y = pred_m)) +
    geom_line(color="white") +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="grey80") +
    ylab("Parasite Prevalence") +
    xlab("Bee community diversity") +
    scale_x_continuous(
        breaks = axis.bee.div,
        labels =  labs.bee.div) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
    #    #theme_classic() +
    theme_dark_black()+
    geom_point(data=spec.all[spec.all$Weights == 1 &
                               spec.all$WeightsPar == 1,],
               aes(y=SiteParasitismRate, x=MeanFloralDiversity, colour="white"))

p6.parasite  <- fit %>%
    spread_draws(b_MeanFloralAbundance_Intercept,
                 b_ParasitePresence_Intercept,
                 b_ParasitePresence_MeanFloralAbundance) %>%
    mutate(MeanFloralAbundance =
               list(seq(min(spec.all$MeanFloralAbundance[
                                           spec.all$WeightsPar == 1]),
                        max(spec.all$MeanFloralAbundance[
                                           spec.all$WeightsPar == 1]),
                        0.1))) %>%
    unnest(MeanFloralAbundance) %>%
    mutate(pred = exp(b_MeanFloralAbundance_Intercept +
                      b_ParasitePresence_Intercept +
                      b_ParasitePresence_MeanFloralAbundance*MeanFloralAbundance)/
               (1+exp(b_MeanFloralAbundance_Intercept +
                      b_ParasitePresence_Intercept +
                      b_ParasitePresence_MeanFloralAbundance*MeanFloralAbundance))) %>%
    group_by(MeanFloralAbundance) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = MeanFloralAbundance, y = pred_m)) +
    geom_line(color="white") +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="grey80") +
    ylab("") +
    xlab("Bee abundance") +
    ylim(0,1)  +
    scale_x_continuous(
        breaks = axis.bee.abund2,
        labels =  labs.bee.abund2) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
            #theme_classic() +
    theme_dark_black()+
    geom_point(data=spec.all[spec.all$Weights == 1 &
                              spec.all$WeightsPar == 1,],
               aes(y=SiteParasitismRate, x=MeanFloralAbundance, colour="white"))

ggsave(p5.parasite, file="figures/parasite_flowerDiv.pdf",
       height=4, width=5)

ggsave(p6.parasite, file="figures/parasite_flowerAbund.pdf",
       height=4, width=5)

parasite.all <- grid.arrange(p5.parasite, p6.parasite, ncol=2)

ggsave(parasite.all, file="figures/flower_parasite.pdf",
       height=4, width=10)


## ***********************************************************************
## bee community- floral div, area
## ***********************************************************************

p1.bee <- fit %>%
    spread_draws(b_PollAbundance_Intercept,
                 b_PollAbundance_Area) %>%
    mutate(Area =
               list(seq(min(spec.all$Area, na.rm=TRUE),
                        max(spec.all$Area,  na.rm=TRUE),
                        0.1))) %>%
    unnest(Area) %>%
    mutate(pred = b_PollAbundance_Intercept +
               b_PollAbundance_Area*Area) %>%
    group_by(Area) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = Area, y = pred_m)) +
    geom_line(color="white") +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="grey80") +
    ylab("Bee abundance") +
    xlab("Meadow area (log)") +
    scale_x_continuous(
        breaks = axis.area,
        labels =  labs.area) +
    scale_y_continuous(
        breaks = axis.bee.abund,
        labels =  labs.bee.abund ) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
        #theme_classic() +
    theme_dark_black()+
    geom_point(data=spec.all[spec.all$Weights == 1,],
               aes(y=PollAbundance, x=Area, color="white"))

p2.bee <- fit %>%
    spread_draws(b_NetBeeDiversity_Intercept, b_MeanFloralDiversity_Intercept,
                 b_NetBeeDiversity_MeanFloralDiversity) %>%
    mutate(MeanFloralDiversity =
               list(seq(min(spec.all$MeanFloralDiversity),
                        max(spec.all$MeanFloralDiversity),
                        0.1))) %>%
    unnest(MeanFloralDiversity) %>%
    mutate(pred = b_NetBeeDiversity_Intercept + b_MeanFloralDiversity_Intercept +
               b_NetBeeDiversity_MeanFloralDiversity*MeanFloralDiversity) %>%
    group_by(MeanFloralDiversity) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = MeanFloralDiversity, y = pred_m)) +
    geom_line(color="white") +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="grey80") +
    ylab("Bee diversity") +
    xlab("Floral diversity") +
    scale_x_continuous(
        breaks = axis.flower.div,
        labels =  labs.flower.div) +
    scale_y_continuous(
        breaks = axis.bee.div,
        labels =  labs.bee.div ) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
        #theme_classic() +
    theme_dark_black()+
    geom_point(data=spec.all[spec.all$Weights == 1,],
               aes(y=PollDiversity, x=MeanFloralDiversity, color="white"))

bee.plots <- grid.arrange(p1.bee, p2.bee, ncol=1)

ggsave(p1.bee, file="figures/bee_area.pdf",
       height=4, width=5)

ggsave(p2.bee, file="figures/bee_floralDiv.pdf",
       height=4, width=5)

ggsave(bee.plots, file="figures/beeComm.pdf",
       height=8, width=4)


## ***********************************************************************
## bee community- area-lat interaction
## ***********************************************************************
## based on "Statistical Rethinking example"

## d  <- spec.uni1
## d <- d %>% as_tibble() %>% select("PollAbundance", "Area", "Lat")
## d <- d[!duplicated(d),]

## d.pivot <- d %>%
##   pivot_wider(
##     names_from = "Lat",
##     values_from = PollAbundance,
##     values_fn = list(PollAbundance = mean)
##   )

## d.pivot <- d.pivot[, c(1, 1 + order(as.numeric(colnames(d.pivot)[-1])))]

## d.pivot  <- data.frame(d.pivot)
## rownames(d.pivot) <- d.pivot$Area
## d.pivot$Area  <- NULL

## df <- scale(d.pivot)
## df[is.na(df)] <- 0

## heatmap(df)

## lats  <- seq(min(spec.uni1$Lat),
##                                    max(spec.uni1$Lat),
##                                    length.out=10)
## ic <-
##     list(Lat = standardize(log(lats)))


## cond.effects <- conditional_effects(fit,
##                                   effects = "Area:Lat",
##                                   int_conditions = ic)

## p.interact  <- plot(cond.effects, plot=FALSE, points=TRUE)

## cols <- viridis(length(ic[[1]]), end=0.9)

## p.interact[[3]] +   ylab("Bee abundance") +
##     xlab("Stand area (log)") +
##     scale_x_continuous(
##         breaks = axis.area,
##         labels =  labs.area) +
##     scale_y_continuous(
##         breaks = axis.bee.abund,
##         labels =  labs.bee.abund ) +
##     theme(axis.title.x = element_text(size=16),
##           axis.title.y = element_text(size=16),
##           text = element_text(size=16)) +
##         #theme_classic() +
    ##theme_dark_black()+
##     scale_colour_discrete(name="Yrs post harvest",
##                           labels=lats,
##                           type=cols) +
##     ## scale_fill_manual(values = alpha(rep("white", length(ic[[1]])),
##     ##                                  1)) +
##     scale_fill_viridis_d(alpha=0.9, option="viridis", end=0.9) +
##     guides(fill = "none")



## ic <-
##   list(Area =  axis.area)

## cond.effects <- conditional_effects(fit,
##                                   effects = "Lat:Area",
##                                   int_conditions = ic)

## p.interact  <- plot(points=T, cond.effects, plot=FALSE)

## cols <- viridis(length(axis.area), end=0.9)

## p.interact[[3]] +   ylab("Bee abundance") +
##     xlab("Years post harvest (log)") +
##     scale_x_continuous(
##         breaks = axis.lat.x,
##         labels =  labs.lat.x) +
##     scale_y_continuous(
##         breaks = axis.bee.abund,
##         labels =  labs.bee.abund ) +
##     theme(axis.title.x = element_text(size=16),
##           axis.title.y = element_text(size=16),
##           text = element_text(size=16)) +
##         #theme_classic() +
   ## theme_dark_black()+
##     scale_colour_discrete(name="Area",
##                           labels=round(exp(labs.area), 0),
##                           type=cols) +
##     scale_fill_viridis_d(alpha=0.01, option="viridis", end=0.9) +
##     guides(fill = "none")


## ***********************************************************************
## bee community- lat
## ***********************************************************************

p2.flower.lat <- fit %>%
    spread_draws(b_MeanFloralDiversity_Intercept,
                 b_MeanFloralDiversity_Lat,
                 ) %>%
    mutate(Lat =
               list(seq(min(spec.all$Lat),
                        max(spec.all$Lat),
                        0.1))) %>%
    unnest(Lat) %>%
    mutate(pred = b_MeanFloralDiversity_Intercept +
               b_MeanFloralDiversity_Lat*Lat) %>%
    group_by(Lat) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = Lat, y = pred_m)) +
    geom_line(color="white") +
    geom_ribbon(aes(ymin =  pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin =  pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="grey80") +
    ylab("Floral Diversity") +
    xlab("Latitude") +
    scale_x_continuous(
        breaks = axis.lat.x,
        labels =  labs.lat.x) +
    scale_y_continuous(
        breaks = axis.flower.div,
        labels =  labs.flower.div) +
    coord_cartesian(ylim = range(axis.flower.div)) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
        #theme_classic() +
    theme_dark_black()+
    geom_point(data=spec.all[spec.all$Weights == 1,],
               aes(y=MeanFloralDiversity, x=Lat, color="white"))

p4.bee.lat <- fit %>%
    spread_draws(b_PollDiversity_Intercept,
                 b_PollDiversity_Lat) %>%
    mutate(Lat =
               list(seq(min(spec.all$Lat),
                        max(spec.all$Lat),
                        0.1))) %>%
    unnest(Lat) %>%
    mutate(pred = b_PollDiversity_Intercept +
               b_PollDiversity_Lat*Lat) %>%
    group_by(Lat) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = Lat, y = pred_m)) +
    geom_line(color="white") +
    geom_ribbon(aes(ymin =  pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin =  pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="grey80") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="grey80") +
    ylab("Bee Diversity") +
    xlab("Latitude") +
    scale_x_continuous(
        breaks = axis.lat.x,
        labels =  labs.lat.x) +
    scale_y_continuous(
        breaks = axis.bee.div,
        labels =  labs.bee.div) +
    coord_cartesian(ylim = range(axis.bee.div)) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
        #theme_classic() +
    theme_dark_black()+
    geom_point(data=spec.all[spec.all$Weights == 1,],
               aes(y=PollDiversity, x=Lat, color="white"))

lat.plots <- grid.arrange(p2.flower.lat,
                          p4.bee.lat,
                          ncol=2)

ggsave(lat.plots, file="figures/lat.pdf",
       height=4, width=8)


p5.parasite.lat <- ggplot(data= spec.all[spec.all$Weights == 1,],
                          aes(x=Lat, y=SiteParasitismRate)) +
    geom_point(color="white")+
    geom_smooth(method=lm) +    theme_dark_black()


ggsave(p5.parasite.lat, file="figures/lat_parasite.pdf",
       height=4, width=5)

