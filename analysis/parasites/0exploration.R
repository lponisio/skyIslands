setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')

## setwd('~/Dropbox (University of Oregon)/skyislands')

setwd("analysis/parasites")
rm(list=ls())

library(ggplot2)
library(lemon)
library(egg)

source("src/misc.R")
source("src/init.R")
source("src/ggplotThemes.R")


site.sum$Year <- as.factor(site.sum$Year)

## ~ Latitute
lat.abund <- ggplot(site.sum, aes(x=Lat, y=PollAbundance, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Bee abundance")

lat.abund.plant <- ggplot(site.sum, aes(x=Lat, y=MeanFloralAbundance, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Flower abundance")

lat.richness <- ggplot(site.sum, aes(x=Lat, y=PollRichness, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Bee richness")

lat.poll.div <- ggplot(site.sum, aes(x=Lat, y=PollDiversity, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Bee diversity")

lat.flower.div <- ggplot(site.sum, aes(x=Lat, y=MeanFloralDiversity, shape=Year, color=Site)) +
    geom_point()+ theme_article() +  theme(legend.position =
                                               "none") + ylab("Flower diversity")

lat.plant.div <- ggplot(site.sum, aes(x=Lat, y=MeanFloweringPlantDiversity, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Flowering plant diversity")

lat.par <- ggplot(site.sum, aes(x=Lat, y=SiteParasitismRate, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Parasite presence")

lat.par.rich <- ggplot(site.sum, aes(x=Lat, y=MeanParasiteRichness, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Parasite richness")

lat <- grid_arrange_shared_legend(lat.abund, lat.abund.plant,
                                  lat.richness, lat.poll.div, lat.flower.div,
                                  lat.plant.div, lat.par, lat.par.rich,
                                  ncol = 2, nrow = 4, position='top')

ggsave(lat, file="../../../skyIslands_saved/figures/Latitute.pdf",
       height=11, width=8.5)


## ~ area
area.abund <- ggplot(site.sum, aes(x=Area, y=PollAbundance, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Bee abundance")

area.abund.plant <- ggplot(site.sum, aes(x=Area, y=MeanFloralAbundance, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Flower abundance")

area.richness <- ggplot(site.sum, aes(x=Area, y=PollRichness, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Bee richness")

area.poll.div <- ggplot(site.sum, aes(x=Area, y=PollDiversity, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Bee diversity")

area.flower.div <- ggplot(site.sum, aes(x=Area, y=MeanFloralDiversity, shape=Year, color=Site)) +
    geom_point()+ theme_article() +  theme(legend.position =
                                               "none") + ylab("Flower diversity")

area.plant.div <- ggplot(site.sum, aes(x=Area, y=MeanFloweringPlantDiversity, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Flowering plant diversity")

area.par <- ggplot(site.sum, aes(x=Area, y=SiteParasitismRate, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Parasite presence")

area.par.rich <- ggplot(site.sum, aes(x=Area, y=MeanParasiteRichness, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Parasite richness")

area <- grid_arrange_shared_legend(area.abund, area.abund.plant,
                                  area.richness, area.poll.div, area.flower.div,
                                  area.plant.div, area.par, area.par.rich,
                                  ncol = 2, nrow = 4, position='top')

ggsave(area, file="../../../skyIslands_saved/figures/area.pdf",
       height=11, width=8.5)


## ~ doy
doy.poll.abund <- ggplot(site.sum, aes(x=SRDoy, y=PollAbundance, color=Site, shape=Year)) +
    geom_point() +
    theme_article() + ylab("Bee abundance") +
    theme(legend.position =  "none")

doy.plant.abund <- ggplot(site.sum, aes(x=SRDoy, y=MeanFloralAbundance, color=Site, shape=Year)) +
    geom_point() +
    theme_article() + ylab("Flower abundance") +
    theme(legend.position =  "none")

doy.poll.div <- ggplot(site.sum, aes(x=SRDoy, y=PollDiversity, color=Site, shape=Year)) +
    geom_point() +
    theme_article()+ ylab("Bee Diversity") +
    theme(legend.position =  "none")

doy.plant.div <- ggplot(site.sum, aes(x=SRDoy, y=MeanFloralDiversity, color=Site, shape=Year)) +
    geom_point() +
    theme_article()+ ylab("Floral Diversity") +
    theme(legend.position =  "none")

doy.par <- ggplot(site.sum, aes(x=SRDoy, y=SiteParasitismRate, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Parasite presence")

doy.par.rich <- ggplot(site.sum, aes(x=SRDoy, y=MeanParasiteRichness, shape=Year, color=Site)) +
    geom_point() + theme_article() + theme(legend.position =
                                               "none") + ylab("Parasite richness")

doy <- grid_arrange_shared_legend(doy.poll.abund, doy.plant.abund,
                                  doy.poll.div, doy.plant.div,
                                  doy.par,doy.par.rich,
                                  ncol = 2, nrow = 3, position='top')

ggsave(doy, file="../../../skyIslands_saved/figures/Doy.pdf",
        height=11, width=8.5)

## max abundance was at SC in 2012, there were a crazy number of tiny
## bees that year
site.sum[site.sum$PollAbundance == max(site.sum$PollAbundance),]



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
