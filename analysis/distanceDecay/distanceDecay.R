## setwd("~/Dropbox/skyIslands/")

rm(list=ls())
wdpath <- 'C:/Users/rah10/Dropbox (University of Oregon)/skyIslands'
setwd(wdpath)
load("data/spec_RBCL_16s.Rdata")

##mantel tests
##https://jkzorz.github.io/2019/07/08/mantel-test.html

##Species abundance dissimilarity matrix: created using a distance measure, i.e. Bray-curtis dissimilarity. This is the same type of dissimilarity matrix used when conducting an ANOSIM test or when making an NMDS plot
#Geographic distance matrix: the physical distance between sites (i.e. Haversine distance)

library(vegan)
library(geosphere)
library(betapart)
library(tidyverse)


meta_cols <- c('UniqueID', 'Family', 'Genus', 'Species', 'GenusSpecies', 'Site', 'Lat', 'Long')

spec16s <- spec.net %>%
  filter(Apidae == 1) %>%
  select(all_of(meta_cols), starts_with('16s')) %>%
  na.omit()

##Species abundance dissimilarity matrix: created using a distance measure, i.e. Bray-curtis dissimilarity. This is the same type of dissimilarity matrix used when conducting an ANOSIM test or when making an NMDS plot
#Geographic distance matrix: the physical distance between sites (i.e. Haversine distance)



## This function takes as input a data frame and a bee genus and filters the data frame to
## just bees of the input genus, computes bray curtis dissimilarity of their 16s communities,
## computes haversine distance matrix of sample sites, performs a mantel test to determine
## correlation between community dissimilarity and distance, then plots the distance decay curves


genusspecies.decay.model <- function(data, type, which){
  #bray curtis dissimilarity matrix of 16s
  
  
  if(type == 'Genus'){
    abund <- data %>%
      filter(Genus == which) %>%
      select(UniqueID, starts_with('16s')) %>%
      select(-UniqueID)
    
    #distance matrix of sites
    geo <- data %>%
      filter(Genus == which) %>%
      select(UniqueID, Long, Lat) %>%
      select(-UniqueID) %>%
      mutate()
  } else if (type == 'GenusSpecies'){
    
    abund <- data %>%
      filter(GenusSpecies == which) %>%
      select(UniqueID, starts_with('16s')) %>%
      select(-UniqueID)
    #distance matrix of sites
    geo <- data %>%
      filter(GenusSpecies == which) %>%
      select(UniqueID, Long, Lat) %>%
      select(-UniqueID) %>%
      mutate()
  }
  
  
  
  
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

bombus_model <- genusspecies.decay.model(spec16s, 'Bombus', type='Genus')
melissodes_model <- genusspecies.decay.model(spec16s, 'Melissodes', type='Genus')

## TODO clean up script
## TODO custom fxn to enable changing axis labs

plot.decay(bombus_model, 
           col='navy', 
           bg=alpha('navy', 0.01), lwd=10,
           cex=2, remove.dots = FALSE,
           xlab='Distance (km)',
           ylab="Bray-Curtis Dissimilarity") 


plot.decay(melissodes_model, 
           col='gold', 
           bg=alpha('gold', 0.1), 
           lwd=10,
           cex=2,
           xlab='Distance (km)',
           ylab="Bray-Curtis Dissimilarity", add=TRUE, remove.dots = FALSE)



