rm(list=ls())

wdpath <- 'C:/Users/rah10/Dropbox (University of Oregon)/PonisioLocal/skyIslands'
setwd(wdpath)

library(tidyverse)
library(igraph)
library(ggplot2)
library(Matrix)

###bee-microbe network attempt using a plant poll network tutorial
### in igraph https://spiesmanecology.com/2017/04/30/network-vis1/

## prep

bee.microbes <- read.csv('spec_RBCL_16s.csv') %>%
  filter(Apidae == 1) #%>%
  # select(UniqueID, starts_with('X16s')) %>%
  # na.omit()

# rownames(bee.microbe.abund) <- bee.microbe.abund$UniqueID
# 
# bee.microbe.abund <- bee.microbe.abund %>%
#   select(-UniqueID) %>%
#   t()
# 
# ## tutorial 
# i_net = graph_from_incidence_matrix(bee.microbe.abund, weight=T)
# V(i_net)$label = NA
# plot(i_net)#,
#      #layout=layout.bipartite)
# 
# 
# ##two mode network
# V(i_net)$color <- c("steel blue", "orange")[V(i_net)$type+1]
# V(i_net)$shape <- c("square", "circle")[V(i_net)$type+1]
# V(i_net)$label <- ""
# V(i_net)$label[V(i_net)$type==F] <- nodes2$media[V(i_net)$type==F]
# V(i_net)$label.cex=.4
# V(i_net)$label.font=2
# 
# #not great
# plot(i_net,
#      vertex.label.color="white",
#      vertex.size=(2-V(i_net)$type)*8,
#      layout=layout.circle)
# 
# # Circle layout -- not great
# l <- layout_in_circle(i_net)
# plot(i_net, layout=l)
# 
# #this is okay but still really busy -- maybe show by genus?
# V(i_net)$size <- degree(i_net)
# V(i_net)$label.cex <- degree(i_net) * 0.2
# plot(i_net, layout = layout_with_graphopt)
# 


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
  plot(i_net, layout = layout.circle)
}

apis_net <- plot_beegenus_network(bee.microbes, 'Apis')

bombus_net <- plot_beegenus_network(bee.microbes, 'Bombus')

anthophora_net <- plot_beegenus_network(bee.microbes, 'Anthophora')

megachile_net <- plot_beegenus_network(bee.microbes, 'Megachile')
