## setwd("~/Dropbox/skyIslands/")
setwd('~/Dropbox (University of Oregon)/skyIslands') ## Rebecca wd
rm(list=ls())
setwd("analysis/distanceDecay")

##mantel tests
##https://jkzorz.github.io/2019/07/08/mantel-test.html

##Species abundance dissimilarity matrix: created using a distance measure, i.e. Bray-curtis dissimilarity. This is the same type of dissimilarity matrix used when conducting an ANOSIM test or when making an NMDS plot
#Geographic distance matrix: the physical distance between sites (i.e. Haversine distance)

library(vegan)
library(geosphere)
library(betapart)
library(tidyverse)


wdpath <- 'C:/Users/rah10/Dropbox (University of Oregon)/PonisioLocal/skyIslands'
setwd(wdpath)

meta_cols <- c('UniqueID', 'Family', 'Genus', 'Species', 'Site', 'Lat', 'Long')

spec16s <- read.csv('spec_RBCL_16s.csv') %>%
  filter(Apidae == 1) %>%
  select(all_of(meta_cols), starts_with('X16s')) %>%
  na.omit()

plot.genus.dissim.dist <- function(data, genus){
#bray curtis dissimilarity matrix
abund <- data %>%
  filter(Genus == genus) %>%
  select(UniqueID, starts_with('X16s')) %>%
  select(-UniqueID)

#distance matrix
geo <- data %>%
  filter(Genus == genus) %>%
  select(UniqueID, Long, Lat) %>%
  select(-UniqueID) %>%
  mutate()

#abundance data frame - bray curtis dissimilarity
dist.abund <- vegdist(abund, method = "bray")

#geographic data frame - haversine distance 
d.geo <- distm(geo, fun = distHaversine)
dist.geo <- as.dist(d.geo)

#abundance vs geographic 
abund_geo  = mantel(dist.abund, dist.geo, method = "spearman", permutations = 9999, na.rm = TRUE)
print(abund_geo)

aa = as.vector(dist.abund)
gg = as.vector(dist.geo)

#new data frame with vectorized distance matrices
mat = data.frame(aa,gg)

#abundance vs geographic distance
# #mm = ggplot(mat, aes(y = aa, x = gg/1000)) + 
#   geom_point(size = 3, alpha = 0.5) + 
#   geom_smooth(method = "lm", colour = "red", alpha = 0.2) + 
#   labs(x = "Physical separation (km)", y = "Bray-Curtis Dissimilarity") + 
#   theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
#          axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
#          axis.title= element_text(face = "bold", size = 14, colour = "black"), 
#          panel.background = element_blank(), 
#          panel.border = element_rect(fill = NA, colour = "black"))
# mm

dist_decay_model <- betapart::decay.model(dist.abund, dist.geo, y.type='dissim')
dist_decay_plot <- plot.decay(dist_decay_model)
dist_decay_plot

}
## by genus plots

plot.genus.dissim.dist(spec16s, 'Apis')
plot.genus.dissim.dist(spec16s, 'Bombus')
plot.genus.dissim.dist(spec16s, 'Anthophora')
plot.genus.dissim.dist(spec16s, 'Megachile')



##################### 1/11/23 distance decay models using betapart decay.model
## https://rdrr.io/cran/betapart/man/decay.model.html

library(betapart)
library(vegan)

dist_decay_model <- betapart::decay.model(dist.abund, dist.geo)


