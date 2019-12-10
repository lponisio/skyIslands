## require(devtools)
## install_github("frankkramer-lab/mully")
setwd('~/Dropbox/skyIslands')
rm(list=ls())
setwd('analysis/multilayer')

library(mully)
library(multinet)
library(ggplot2)
library(betalink)
library(igraph)

load('../../data/spec.Rdata')
load('../../data/nets.Rdata')

site.order <- c("JC", "SC", "SM", "MM", "HM", "PL", "CH")

years <- sapply(strsplit(names(nets), "[.]"), function(x) x[[2]])
sites <- sapply(strsplit(names(nets), "[.]"), function(x) x[[1]])

nets.graph <- prepare_networks(nets, directed=FALSE)
nets.graph <- nets.graph

same.year <- split(nets.graph, years)

## **************************************************************
## 2012 multi layer network
## **************************************************************
net.2012 <- ml_empty()
add_igraph_layer_ml(net.2012, same.year$'2012'$JC.2012.1, "layer1")
add_igraph_layer_ml(net.2012, same.year$'2012'$SC.2012.1, "layer2")
add_igraph_layer_ml(net.2012, same.year$'2012'$MM.2012.1, "layer3")
add_igraph_layer_ml(net.2012, same.year$'2012'$PL.2012.1, "layer4")
add_igraph_layer_ml(net.2012, same.year$'2012'$CH.2012.1, "layer5")

l1.to.l2 <- V(same.year$'2012'$JC.2012.1)[V(same.year$'2012'$JC.2012.1) %in%
                              V(same.year$'2012'$SC.2012.1)]
l1.to.l2 <- names(l1.to.l2)

l2.to.l3 <- V(same.year$'2012'$SC.2012.1)[V(same.year$'2012'$SC.2012.1) %in%
                              V(same.year$'2012'$MM.2012.1)]
l2.to.l3 <- names(l2.to.l3)

l3.to.l4 <- V(same.year$'2012'$MM.2012.1)[V(same.year$'2012'$MM.2012.1) %in%
                              V(same.year$'2012'$PL.2012.1)]
l3.to.l4 <- names(l3.to.l4)

l4.to.l5 <- V(same.year$'2012'$PL.2012.1)[V(same.year$'2012'$PL.2012.1) %in%
                              V(same.year$'2012'$CH.2012.1)]
l4.to.l5 <- names(l4.to.l5)

inter_layer_edges12 <- data.frame(
    actor1 = l1.to.l2,
    layer1 = rep("layer1", length(l1.to.l2)),
    actor2 =  l1.to.l2,
    layer2 = rep("layer2", length(l1.to.l2)))

add_edges_ml(net.2012, inter_layer_edges12)

inter_layer_edges23 <- data.frame(
    actor1 = l2.to.l3,
    layer1 = rep("layer2", length(l2.to.l3)),
    actor2 =  l2.to.l3,
    layer2 = rep("layer3", length(l2.to.l3)))

add_edges_ml(net.2012, inter_layer_edges23)

## inter_layer_edges34 <- data.frame(
##     actor1 = l3.to.l4,
##     layer1 = rep("layer3", length(l3.to.l4)),
##     actor2 =  l3.to.l4,
##     layer2 = rep("layer4", length(l3.to.l4)))

## add_edges_ml(net.2012, inter_layer_edges34)

## inter_layer_edges45 <- data.frame(
##     actor1 = l4.to.l5,
##     layer1 = rep("layer4", length(l4.to.l5)),
##     actor2 =  l4.to.l5,
##     layer2 = rep("layer5", length(l4.to.l5)))

## add_edges_ml(net.2012, inter_layer_edges45)


l <- layout_multiforce_ml(net.2012, w_inter = 1, gravity = 1)
bk <- par("mar")
par(mar=c(0,0,0,0))
plot(net.2012, layers=c("layer1", "layer2", "layer3", "layer4",
                        "layer5"),
                        layout = l, vertex.labels = "",
     legend.x = "bottomright", legend.inset = c(.05, .05))
par(mar=bk)
