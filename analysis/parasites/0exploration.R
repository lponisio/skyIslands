## setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')

setwd('~/Dropbox (University of Oregon)/skyislands')

setwd("analysis/parasites")
rm(list=ls())

library(ggplot2)
library(lemon)
library(egg)

source("src/misc.R")
source("src/init.R")
source("src/ggplotThemes.R")


site.sum$Year <- as.factor(site.sum$Year)
## site.sum <- site.sum[site.sum$Year != 2022,]

## ********************************************
## Bombus abundance relationships
## ********************************************

hb.bombus <- ggplot(site.sum, aes(x=BombusAbundance,
                                  y=HBAbundance, shape=Year, color=Site)) +
    geom_point() + theme_article() +
    ylab("Apis abundance") +
    xlab("Bombus abundance") +
    theme(legend.position = "none")
bombus.nonbombusHB <- ggplot(site.sum, aes(x=BombusAbundance,
                                           y=NonBombusHBAbundance,
                                           shape=Year, color=Site)) +
    geom_point() + theme_article() +
    xlab("Bombus abundance") +
    ylab("Non-Apis abundance") +
    theme(legend.position =  "none")
bombus.polldiv <- ggplot(site.sum, aes(x=BombusAbundance,
                                       y=PollDiversity,
                                       shape=Year, color=Site)) +
    geom_point() + theme_article() +
    xlab("Bombus abundance") +
    ylab("Pollinator diversity") +
    theme(legend.position =  "none")
bombus.comps <- grid_arrange_shared_legend(hb.bombus, bombus.nonbombusHB,
                                  bombus.polldiv,
                                  ncol = 1, nrow = 3, position='top',
                                  plot=FALSE)
ggsave(bombus.comps, file="../../../skyIslands_saved/figures/bombus_comps.pdf",
       height=8, width=4)

## ********************************************
## Apis abundance relationships
## ********************************************

HB.nonHB <- ggplot(site.sum, aes(x=HBAbundance,
                                           y=NonBombusHBAbundance,
                                           shape=Year, color=Site)) +
    geom_point() + theme_article() +
    xlab("Apis abundance") +
    ylab("Non-Apis/Bombus abundance") +
    theme(legend.position =  "none")
HB.polldiv <- ggplot(site.sum, aes(x=HBAbundance,
                                       y=PollDiversity,
                                       shape=Year, color=Site)) +
    geom_point() + theme_article() +
    xlab("Apis abundance") +
    ylab("Pollinator diversity") +
    theme(legend.position =  "none")
apis.comps <- grid_arrange_shared_legend(HB.nonHB,
                                  HB.polldiv,
                                  ncol = 1, nrow = 2,
                                  position='top',
                                   plot=FALSE)
ggsave(apis.comps, file="../../../skyIslands_saved/figures/apis_comps.pdf",
       height=8, width=4)

## ********************************************
## Insect Abundance ~ Latitude
## ********************************************

lat.abund <- ggplot(site.sum, aes(x=Lat, y=PollAbundance, shape=Year, color=Site)) +
    geom_point() + theme_article() +
    theme(legend.position = "none") +
    ylab("Bee abundance") +
    xlab("Latitude")
lat.abund.hb  <- ggplot(site.sum, aes(x=Lat, y=HBAbundance, shape=Year, color=Site)) +
    geom_point() + theme_article() +
    theme(legend.position = "none") +
    ylab("Apis abundance") +
    xlab("Latitude")
lat.abund.bombus  <- ggplot(site.sum,
                            aes(x=Lat, y=BombusAbundance, shape=Year, color=Site)) +
    geom_point() + theme_article() +
    theme(legend.position = "none") +
    ylab("Bombus abundance") +
    xlab("Latitude")
lat.abund.other  <- ggplot(site.sum,
                           aes(x=Lat, y=NonBombusHBAbundance, shape=Year, color=Site)) +
    geom_point() + theme_article() +
    theme(legend.position = "none") +
    ylab("Non-Apis/Bombus abundance") +
    xlab("Latitude")
lat.abund.all <- grid_arrange_shared_legend(lat.abund, lat.abund.hb,
                                  lat.abund.bombus, lat.abund.other,
                                  ncol = 2, nrow = 2, position='top',
                                  plot=FALSE)
ggsave(lat.abund.all, file="../../../skyIslands_saved/figures/Latitute_abundance.pdf",
       height=9, width=8.5)

## ********************************************
## Insect Diversity ~ Latitude
## ********************************************

lat.div <- ggplot(site.sum, aes(x=Lat, y=PollDiversity, shape=Year, color=Site)) +
    geom_point() + theme_article() +
    theme(legend.position = "none") +
    ylab("Bee diversity") +
    xlab("Latitude")
lat.rich  <- ggplot(site.sum, aes(x=Lat, y=PollRichness, shape=Year, color=Site)) +
    geom_point() + theme_article() +
    theme(legend.position = "none") +
    ylab("Pollinator richness") +
    xlab("Latitude")
lat.div.bombus  <- ggplot(site.sum,
                            aes(x=Lat, y=BombusDiversity, shape=Year, color=Site)) +
    geom_point() + theme_article() +
    theme(legend.position = "none") +
    ylab("Bombus diversity") +
    xlab("Latitude")
lat.rich.bombus  <- ggplot(site.sum,
                            aes(x=Lat, y=BombusRichness, shape=Year, color=Site)) +
    geom_point() + theme_article() +
    theme(legend.position = "none") +
    ylab("Bombus richness") +
    xlab("Latitude")
lat.div.all <- grid_arrange_shared_legend(lat.div, lat.rich,
                                  lat.div.bombus, lat.rich.bombus,
                                  ncol = 2, nrow = 2, position='top',
                                  plot=FALSE)
ggsave(lat.div.all, file="../../../skyIslands_saved/figures/Latitute_diversity.pdf",
       height=9, width=8.5)

## ********************************************
## Plants ~ Latitude
## ********************************************

lat.abund.plant <- ggplot(site.sum,
                          aes(x=Lat, y=MeanFloralAbundance,
                              shape=Year,
                              color=Site)) +
    geom_point() + theme_article() +
    theme(legend.position = "none") +
    ylab("Flower abundance")
lat.plant.div <- ggplot(site.sum,
                        aes(x=Lat,
                            y=MeanFloweringPlantDiversity,
                            shape=Year,
                            color=Site)) +
    geom_point() + theme_article() +
    theme(legend.position = "none") +
    ylab("Flowering plant diversity")
lat.flower.div <- ggplot(site.sum,
                         aes(x=Lat, y=MeanFloralDiversity,
                             shape=Year,
                             color=Site)) +
    geom_point()+ theme_article() +
    theme(legend.position = "none") +
    ylab("Flower diversity")
lat.visit.div <- ggplot(site.sum,
                        aes(x=Lat, y=VisitedFloralDiversity,
                            shape=Year,
                            color=Site)) +
    geom_point()+ theme_article() +
    theme(legend.position = "none") +
    ylab("Visited flower diversity")
lat.flower.rich <- ggplot(site.sum,
                          aes(x=Lat, y=MeanFloralRichness,
                              shape=Year,
                              color=Site)) +
    geom_point()+ theme_article() +
    theme(legend.position = "none") +
    ylab("Flower richness")
lat.visit.rich <- ggplot(site.sum,
                         aes(x=Lat, y=VisitedFloralRichness,
                             shape=Year,
                             color=Site)) +
    geom_point()+ theme_article() +
    theme(legend.position = "none") +
    ylab("Visited flower richness")
lat.plant <- grid_arrange_shared_legend(lat.abund.plant, lat.plant.div,
                                        lat.flower.div, lat.visit.div,
                                        lat.flower.rich, lat.visit.rich,
                                        ncol = 2, nrow = 3,
                                        position='top',
                                        plot=FALSE)

ggsave(lat.plant, file="../../../skyIslands_saved/figures/Latitute_plants.pdf",
       height=11, width=8.5)

## ********************************************
## Parasites ~ Latitude
## ********************************************

lat.par <- ggplot(site.sum,
                  aes(x=Lat, y=SiteParasitismRate, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position ="none") +
    ylab("Parasite presence")
lat.par.rich <- ggplot(site.sum,
                       aes(x=Lat, y=MeanParasiteRichness, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position = "none") +
    ylab("Parasite richness")
lat.par.all <- grid_arrange_shared_legend(lat.par, lat.par.rich,
                                          ncol = 1, nrow = 2, position='top',
                                          plot=FALSE)
ggsave(lat.par.all, file="../../../skyIslands_saved/figures/Latitute_parasites.pdf",
       height=11, width=8.5)

## ********************************************
## ~ Area
## ********************************************

area.abund <- ggplot(site.sum,
                     aes(x=Area, y=PollAbundance, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Bee abundance")
area.abund.plant <- ggplot(site.sum,
                           aes(x=Area, y=MeanFloralAbundance, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Flower abundance")
area.richness <- ggplot(site.sum,
                        aes(x=Area, y=PollRichness, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Bee richness")
area.poll.div <- ggplot(site.sum,
                        aes(x=Area, y=PollDiversity, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Bee diversity")
area.flower.div <- ggplot(site.sum,
                          aes(x=Area, y=MeanFloralDiversity, shape=Year, color=Site)) +
    geom_point()+ theme_article() +  theme(legend.position =
                                               "none") + ylab("Flower diversity")
area.plant.div <- ggplot(site.sum,
                         aes(x=Area, y=MeanFloweringPlantDiversity, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Flowering plant diversity")

area.par <- ggplot(site.sum,
                   aes(x=Area, y=SiteParasitismRate, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Parasite presence")

area.par.rich <- ggplot(site.sum,
                        aes(x=Area, y=MeanParasiteRichness, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Parasite richness")
area <- grid_arrange_shared_legend(area.abund, area.abund.plant,
                                  area.richness, area.poll.div, area.flower.div,
                                  area.plant.div, area.par, area.par.rich,
                                  ncol = 2, nrow = 4, position='top',
                                  plot=FALSE)
ggsave(area, file="../../../skyIslands_saved/figures/area.pdf",
       height=11, width=8.5)

## ********************************************
## ~ doy
## ********************************************

doy.poll.abund <- ggplot(site.sum,
                         aes(x=SRDoy, y=PollAbundance, color=Site, shape=Year)) +
    geom_point() +
    theme_article() + ylab("Bee abundance") +
    theme(legend.position =  "none")

doy.HB.abund <- ggplot(site.sum,
                       aes(x=SRDoy, y=HBAbundance, color=Site, shape=Year)) +
    geom_point() +
    theme_article() + ylab("Apis abundance") +
    theme(legend.position =  "none")


doy.bombus.abund <- ggplot(site.sum,
                           aes(x=SRDoy, y=BombusAbundance, color=Site, shape=Year)) +
    geom_point() +
    theme_article() + ylab("Bombus abundance") +
    theme(legend.position =  "none")


doy.plant.abund <- ggplot(site.sum,
                          aes(x=SRDoy, y=MeanFloralAbundance, color=Site, shape=Year)) +
    geom_point() +
    theme_article() + ylab("Flower abundance") +
    theme(legend.position =  "none")

doy.poll.div <- ggplot(site.sum,
                       aes(x=SRDoy, y=PollDiversity, color=Site, shape=Year)) +
    geom_point() +
    theme_article()+ ylab("Bee Diversity") +
    theme(legend.position =  "none")

doy.bombus.div <- ggplot(site.sum,
                       aes(x=SRDoy, y=BombusDiversity, color=Site, shape=Year)) +
    geom_point() +
    theme_article()+ ylab("Bombus Diversity") +
    theme(legend.position =  "none")

doy.plant.div <- ggplot(site.sum,
                        aes(x=SRDoy, y=MeanFloralDiversity, color=Site, shape=Year)) +
    geom_point() +
    theme_article()+ ylab("Floral Diversity") +
    theme(legend.position =  "none")

doy.par <- ggplot(site.sum,
                  aes(x=SRDoy, y=SiteParasitismRate, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Parasite presence")

doy.par.rich <- ggplot(site.sum,
                       aes(x=SRDoy,  y=MeanParasiteRichness, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Parasite richness")
doy <- grid_arrange_shared_legend(doy.poll.abund,
                                  doy.HB.abund,
                                  doy.bombus.abund,
                                  doy.plant.abund,
                                  doy.poll.div,
                                  doy.bombus.div,
                                  doy.plant.div,
                                  doy.par,doy.par.rich,
                                  ncol = 3, nrow = 3, position='top',
                                  plot=FALSE)

ggsave(doy, file="../../../skyIslands_saved/figures/Doy.pdf",
        height=11, width=8.5)

## max abundance was at SC in 2012, there were a crazy number of tiny
## bees that year
site.sum[site.sum$PollAbundance == max(site.sum$PollAbundance),]

## ********************************************
## Parasite rate barchart
## ********************************************

par.counts <- table(spec$ParasiteRichness[spec$Apidae == 1])

par.counts <- par.counts/sum(par.counts)

par.counts <- as.data.frame(par.counts)

p <- ggplot(data=par.counts, aes(x=Var1, y=Freq)) +
    geom_bar(stat="identity", fill="darkolivegreen") +
    theme_dark_black() + coord_flip() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold")) +
              xlab("Number of parasites/individual") +
              ylab("Proportion of screened bees")

ggsave(p, file="../../../skyIslands_saved/figures/summary_parasites.pdf",
        height=4, width=5)

## ********************************************
## Site characteristics
## ********************************************

elev.lat <- ggplot(site.sum,
                   aes(x=Lat,  y=Elev, color=Site)) +
    geom_point() + theme_article() +
    theme(legend.position = "none") +
    ylab("Elevation") +
    xlab("Latitude")


area.lat <- ggplot(site.sum,
                   aes(x=Lat,  y=Area, color=Site)) +
    geom_point() + theme_article() +
    theme(legend.position = "none") +
    ylab("Area") +
    xlab("Latitude")

site <- grid_arrange_shared_legend(elev.lat,
                                   area.lat,
                                  ncol = 2, nrow = 1, position='top',
                                  plot=FALSE)


ggsave(site, file="../../../skyIslands_saved/figures/site_char.pdf",
        height=3, width=6)
