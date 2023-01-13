rm(list=ls())

wdpath <- 'C:/Users/rah10/Dropbox (University of Oregon)/PonisioLocal/skyIslands'
setwd(wdpath)

library(tidyverse)
library(igraph)
library(ggplot2)

###bee-microbe network attempt using a plant poll network tutorial
### in igraph https://spiesmanecology.com/2017/04/30/network-vis1/

## prep

bee.microbe.abund <- read.csv('spec_RBCL_16s.csv') %>%
  filter(Apidae == 1) %>%
  select(UniqueID, starts_with('X16s')) %>%
  na.omit()

rownames(bee.microbe.abund) <- bee.microbe.abund$UniqueID

bee.microbe.abund <- bee.microbe.abund %>%
  select(-UniqueID) %>%
  t()

## tutorial 
i_net = graph_from_incidence_matrix(bee.microbe.abund, weight=T)
V(i_net)$label = NA
plot(i_net)#,
     #layout=layout.bipartite)


