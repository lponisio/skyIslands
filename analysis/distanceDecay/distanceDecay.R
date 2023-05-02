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


genus.decay.model <- function(data, genus){
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
dist.geo <- as.dist(d.geo)/1000

#abundance vs geographic mantel test
abund_geo  = mantel(dist.abund, dist.geo, method = "spearman", permutations = 9999, na.rm = TRUE)
print(abund_geo)

dist_decay_model <- betapart::decay.model(dist.abund,
                                          dist.geo,
                                          y.type='dissim',
                                          model.type = 'exp',
                                          perm=100)
# dist_decay_plot <- plot.decay(dist_decay_model,
#                               main=genus)
# dist_decay_plot
dist_decay_model

}
## by genus plots

apis_model <- genus.decay.model(spec16s, 'Apis')

ci95 <- predict(apis_model$model, interval = "confidence", level = 0.95)

bombus_model <- genus.decay.model(spec16s, 'Bombus')
anthophora_model <- genus.decay.model(spec16s, 'Anthophora')
megachile_model <- genus.decay.model(spec16s, 'Megachile')

plot.decay(bombus_model, 
           col='#ED7953', 
           bg=alpha('#ED7953', 0.1), 
           pch = 22, lwd=10,
           cex=2) 

plot.decay(apis_model, 
           col='#9C179E', 
           bg=alpha('#9C179E', 0.1), 
           pch = 21, lwd=10,
           cex=2,
           xlab='Distance (km)',
           ylab="Bray-Curtis Dissimilarity", add=TRUE) 


plot.decay(megachile_model, 
           col='#F0F921', 
           bg=alpha('#F0F921', 0.1), 
           pch = 23, lwd=10,
           cex=2,
           xlab='Distance (km)',
           ylab="Bray-Curtis Dissimilarity", add=TRUE) 

plot.decay(anthophora_model, 
           col='#0D0887', 
           bg=alpha('#0D0887', 0.1), 
           pch = 24, lwd=10,
           cex=2,
           xlab='Distance (km)',
           ylab="Bray-Curtis Dissimilarity", add=TRUE) 