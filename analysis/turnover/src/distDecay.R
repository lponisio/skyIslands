##mantel tests
##https://jkzorz.github.io/2019/07/08/mantel-test.html

##Species abundance dissimilarity matrix: created using a distance measure, i.e. Bray-curtis dissimilarity. This is the same type of dissimilarity matrix used when conducting an ANOSIM test or when making an NMDS plot
#Geographic distance matrix: the physical distance between sites (i.e. Haversine distance)

library(vegan)
library(geosphere)
library(betapart)
library(tidyverse)

## **********************************************************
## Standardize, center, and transform data
## **********************************************************

## all of the variables that are explanatory variables and thus need
## to be centered

variables.to.log <- c("rare.degree",
                      "BeeAbundance"
)

vars_yearsr <- c("MeanFloralAbundance",
                 "MeanFloralDiversity",
                 "SRDoy",
                 "BeeAbundance",
                 "BeeDiversity",
                 "VisitedFloralDiversity"
)

vars_yearsrsp <- "rare.degree"
vars_sp <- "MeanITD.x"
vars_site <- "Lat"

source("../microbiome/src/misc_microbe.R")
source("../microbiome/src/standardize_weights_parasites.R")
source("../microbiome/src/makeMultiLevelData.R")
source("../microbiome/src/init_microbe.R")


meta_cols <- c('UniqueID', 'Family', 'Genus', 'Species', 'GenusSpecies', 'Site', 'Lat', 'Long', 'WeightsObligateMicrobe', 'WeightsTransientMicrobe')

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


genusspecies.decay.model <- function(data, type, which, model.type){
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
  abund_geo  = mantel(dist.abund, dist.geo, method = "spearman", permutations = 999, na.rm = TRUE)
  print(abund_geo)
  
  dist_decay_model <- betapart::decay.model(dist.abund,
                                            dist.geo,
                                            y.type='dissim',
                                            model.type = model.type,
                                            perm=999)
  # dist_decay_plot <- plot.decay(dist_decay_model,
  #                               main=genus)
  # dist_decay_plot
  dist_decay_model
  
}


microbe.type.decay.model <- function(data, type, model.type){
  #bray curtis dissimilarity matrix of 16s
  
  
  if(type == 'Obligate'){
    abund <- data %>%
      filter(WeightsObligateMicrobe == 1) %>%
      select(UniqueID, starts_with('16s')) %>%
      select(-UniqueID)
    
    #distance matrix of sites
    geo <- data %>%
      filter(WeightsObligateMicrobe == 1) %>%
      select(UniqueID, Long, Lat) %>%
      select(-UniqueID) %>%
      mutate()
  } else if (type == 'Facultative'){
    abund <- data %>%
      filter(WeightsTransientMicrobe == 1) %>%
      select(UniqueID, starts_with('16s')) %>%
      select(-UniqueID)
    
    #distance matrix of sites
    geo <- data %>%
      filter(WeightsTransientMicrobe == 1) %>%
      select(UniqueID, Long, Lat) %>%
      select(-UniqueID) %>%
      mutate()
  }
  
  
  
  
  #abundance data frame - bray curtis dissimilarity
  dist.abund <- vegdist(abund, method = "bray")
  
  #geographic data frame - haversine distance in m (takes a df with lat and long and calculates dist)
  d.geo <- distm(geo, fun = distHaversine)
  dist.geo <- as.dist(d.geo)/1000
  
  #browser()
  
  #abundance vs geographic mantel test
  abund_geo  = mantel(dist.abund, dist.geo, method = "spearman", permutations = 999, na.rm = TRUE)
  print(abund_geo)
  
  dist_decay_model <- betapart::decay.model(dist.abund,
                                            dist.geo,
                                            y.type='dissim',
                                            model.type = model.type,
                                            perm=999)
  # dist_decay_plot <- plot.decay(dist_decay_model,
  #                               main=genus)
  # dist_decay_plot
  dist_decay_model
  
}


plot_decay_ggplot_single <- function(x,
                                     xlab,
                                     ylab,
                                     mod1color='navy',
                                     col = "black",
                                     lty = "solid",
                                     lwd = 1.5,
                                     cex = 1) {
  #browser()
  # Extract data and fitted values
  data <- data.frame(x$data.x, x$data.y)
  model <- x$model
  fitted_values <- data.frame(fitted(model)[order(data$x.data.x)])
  sorted_data <- data[order(data$x.data.x), ]
  #browser()
  # Create the ggplot object
  p <- ggplot(data, aes(x = x.data.x, y = x.data.y)) +
    geom_jitter(fill = mod1color, alpha = 0.1 , color="black", pch=21, cex=3) +
    geom_line(aes(x = sorted_data$x.data.x, y = (1 - fitted_values$fitted.model..order.data.x.data.x..)),
              color = 'black', linewidth=2.5,) +
    geom_line(aes(x = sorted_data$x.data.x, y = (1 - fitted_values$fitted.model..order.data.x.data.x..)), color = mod1color, linewidth=2) +
    labs(x=xlab, y=ylab) +
    theme_classic() +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16))
  
  return(p)
}


plot_decay_ggplot_combined <- function(x,
                                       z,
                                       xlab,
                                       ylab,
                                       mod1color='navy',
                                       mod2color='gold',
                                       alpha1 = 0.1,
                                       alpha2 = 0.5,
                                       col = "black",
                                       lty1,
                                       lty2,
                                       lwd = 1.5,
                                       cex = 1) {
  
  # Extract data and fitted values
  data1 <- data.frame(x$data.x, x$data.y)
  model1 <- x$model
  fitted_values1 <- data.frame(fitted(model1)[order(data1$x.data.x)])
  sorted_data1 <- data1[order(data1$x.data.x), ]
  
  # Extract data and fitted values
  data2 <- data.frame(z$data.x, z$data.y)
  model2 <- z$model
  fitted_values2 <- data.frame(fitted(model2)[order(data2$z.data.x)])
  sorted_data2 <- data2[order(data2$z.data.x), ]
  
  #browser()
  
  # Create the ggplot object
  p <- ggplot(data1, aes(x = x.data.x, y = x.data.y)) +
    geom_point(data=data1, aes(x = x.data.x, y = x.data.y), fill = mod1color, alpha = alpha1 , color="black", pch=21, cex=3, position = position_jitter(w = 0.5, h = 0)) +
    geom_point(data=data2, aes(x = z.data.x, y = z.data.y), fill = mod2color, alpha = alpha2 , color="black", pch=21, cex=3, position = position_jitter(w = 0.5, h = 0)) +
    geom_line(aes(x = sorted_data1$x.data.x, y = (1 - fitted_values1$fitted.model1..order.data1.x.data.x..)),
              color = 'black', linewidth=2.5,) +
    geom_line(aes(x = sorted_data1$x.data.x, y = (1 - fitted_values1$fitted.model1..order.data1.x.data.x..)), color = mod1color, linewidth=2, linetype=lty1) +
    geom_line(data=data2, aes(x = sorted_data2$z.data.x, y = (1 - fitted_values2$fitted.model2..order.data2.z.data.x..)),
              color = 'black', linewidth=2.5,) +
    geom_line(data=data2, aes(x = sorted_data2$z.data.x, y = (1 - fitted_values2$fitted.model2..order.data2.z.data.x..)), color = mod2color, linewidth=2, linetype=lty2) +
    labs(x=xlab, y=ylab) +
    theme_classic() +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16))
  
  return(p)
}
