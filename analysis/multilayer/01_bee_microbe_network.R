rm(list=ls())

wdpath <- 'C:/Users/rah10/Dropbox (University of Oregon)/PonisioLocal/skyIslands'
setwd(wdpath)

library(tidyverse)
library(igraph)
library(ggplot2)
library(Matrix)
library(ggpubr)

###bee-microbe network attempt using a plant poll network tutorial
### in igraph https://spiesmanecology.com/2017/04/30/network-vis1/

## prep

bee.microbes <- read.csv('spec_RBCL_16s.csv') %>%
  filter(Apidae == 1) #%>%

#function to plot bee-microbe network by bee genus

plot_beegenus_network <- function(data, beegenus){
  bee.microbe.abund <- data %>%
    filter(Apidae == 1) %>%
    filter(Genus == beegenus) %>%
    select(UniqueID, starts_with('X16s')) %>%
    na.omit()
  
  rownames(bee.microbe.abund) <- bee.microbe.abund$UniqueID
  
  bee.microbe.abund <- bee.microbe.abund %>%
    select(-UniqueID) %>%
    t()
  
  i_net = graph_from_incidence_matrix(bee.microbe.abund, weight=T)

  V(i_net)$label = NA
  
  V(i_net)$color <- c("steel blue", "orange")[V(i_net)$type+1]
  V(i_net)$shape <- c("square", "circle")[V(i_net)$type+1]
  V(i_net)$size <- degree(i_net)
  V(i_net)$label.cex <- degree(i_net) * 0.1
  plot(i_net, layout = layout.circle) + title(main=beegenus)
}

apis_net <- plot_beegenus_network(bee.microbes, 'Apis')

bombus_net <- plot_beegenus_network(bee.microbes, 'Bombus')

anthophora_net <- plot_beegenus_network(bee.microbes, 'Anthophora')

megachile_net <- plot_beegenus_network(bee.microbes, 'Megachile')


##### trying now to make igraph objects for each site

find_site_network <- function(data, site){
  bee.microbe.abund <- data %>%
    filter(Apidae == 1) %>%
    filter(Site == site) %>%
    select(UniqueID, starts_with('X16s')) %>%
    na.omit()
  
  rownames(bee.microbe.abund) <- bee.microbe.abund$UniqueID
  
  bee.microbe.abund <- bee.microbe.abund %>%
    select(-UniqueID) %>%
    t()
  
  i_net = graph_from_incidence_matrix(bee.microbe.abund, weight=T)
  
  #i_net
  V(i_net)$label = NA

  V(i_net)$color <- c("steel blue", "orange")[V(i_net)$type+1]
  V(i_net)$shape <- c("square", "circle")[V(i_net)$type+1]
  V(i_net)$size <- degree(i_net)
  V(i_net)$label.cex <- degree(i_net) * 0.1
  plot(i_net, layout = layout.circle) + title(main=site)
  return(i_net)
}

site_list <- c('JC', 'SM', 'SC', 'MM', 'HM', 'PL', 'CH')

#network_list <- array() ##make a list of all networks

for (i in site_list){
  this_network <- find_site_network(bee.microbes, i)
}


## put several webs into an array, blowing the dimensions up to the union of species -- not working need to use for loop to generate 7 distinct webs then add them in a list
all_nets <- webs2array(network_list[1],
                       network_list[2],
                       network_list[3])

##use betalinkr_multi once list is fixed





