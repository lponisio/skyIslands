## setwd("~/Dropbox/skyIslands/")
setwd('~/Dropbox (University of Oregon)/skyIslands') ## Rebecca wd
rm(list=ls())
setwd("analysis/distanceDecay")

##mantel tests
##https://jkzorz.github.io/2019/07/08/mantel-test.html

##Species abundance dissimilarity matrix: created using a distance measure, i.e. Bray-curtis dissimilarity. This is the same type of dissimilarity matrix used when conducting an ANOSIM test or when making an NMDS plot
#Geographic distance matrix: the physical distance between sites (i.e. Haversine distance)

library(vegan)
install.packages("geosphere")
library(geosphere)


wdpath <- 'C:/Users/rah10/Dropbox (University of Oregon)/PonisioLocal/skyIslands'
setwd(wdpath)

meta_cols <- c('UniqueID', 'Family', 'Genus', 'Species', 'Site', 'Lat', 'Long')

spec16s <- read.csv('spec_RBCL_16s.csv') %>%
  filter(Apidae == 1) %>%
  select(all_of(meta_cols), starts_with('X16s'))


#bray curtis dissimilarity matrix
bray <- spec16s %>%
  select(UniqueID, starts_with('X16s'))

#distance matrix
dist0 <- spec16s %>%
  select(UniqueID, Lat, Long)

