## setwd("~/Dropbox/skyIslands/")

rm(list=ls())

wdpath <- 'C:/Users/rah10/Dropbox (University of Oregon)/PonisioLocal/skyIslands'
setwd(wdpath)

##mantel tests
##https://jkzorz.github.io/2019/07/08/mantel-test.html

##Species abundance dissimilarity matrix: created using a distance measure, i.e. Bray-curtis dissimilarity. This is the same type of dissimilarity matrix used when conducting an ANOSIM test or when making an NMDS plot
#Geographic distance matrix: the physical distance between sites (i.e. Haversine distance)

library(vegan)
library(geosphere)
library(betapart)
library(tidyverse)


meta_cols <- c('UniqueID', 'Family', 'Genus', 'Species', 'Site', 'Lat', 'Long')

spec16s <- read.csv('spec_RBCL_16s.csv') %>%
  filter(Apidae == 1) %>%
  select(all_of(meta_cols), starts_with('X16s')) %>%
  na.omit()

##Species abundance dissimilarity matrix: created using a distance measure, i.e. Bray-curtis dissimilarity. This is the same type of dissimilarity matrix used when conducting an ANOSIM test or when making an NMDS plot
#Geographic distance matrix: the physical distance between sites (i.e. Haversine distance)



## This function takes as input a data frame and a bee genus and filters the data frame to
## just bees of the input genus, computes bray curtis dissimilarity of their 16s communities,
## computes haversine distance matrix of sample sites, performs a mantel test to determine
## correlation between community dissimilarity and distance, then plots the distance decay curves


plot.genus.dist.decay <- function(data, genus){
#bray curtis dissimilarity matrix of 16s
abund <- data %>%
  filter(Genus == genus) %>%
  select(UniqueID, starts_with('X16s')) %>%
  select(-UniqueID)

#distance matrix of sites
geo <- data %>%
  filter(Genus == genus) %>%
  select(UniqueID, Long, Lat) %>%
  select(-UniqueID) %>%
  mutate()

#abundance data frame - bray curtis dissimilarity
dist.abund <- vegdist(abund, method = "bray")

#geographic data frame - haversine distance in m (takes a df with lat and long and calculates dist)
d.geo <- distm(geo, fun = distHaversine)
dist.geo <- as.dist(d.geo)

#abundance vs geographic mantel test
abund_geo  = mantel(dist.abund, dist.geo, method = "spearman", permutations = 9999, na.rm = TRUE)
print(abund_geo)

dist_decay_model <- betapart::decay.model(dist.abund, dist.geo, y.type='dissim')
dist_decay_plot <- plot.decay(dist_decay_model, main=genus)
dist_decay_plot

}
## by genus plots

plot.genus.dist.decay(spec16s, 'Apis')
plot.genus.dist.decay(spec16s, 'Bombus')
plot.genus.dist.decay(spec16s, 'Anthophora')
plot.genus.dist.decay(spec16s, 'Megachile')





