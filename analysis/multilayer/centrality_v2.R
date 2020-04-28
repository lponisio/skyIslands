setwd('~/Dropbox/skyIslands')
rm(list=ls())
setwd('analysis/multilayer')
library(igraph)

net.type <- "Yr"
species <- c("Plant", "Pollinator")
## species <- c("Pollinator", "Parasite")

source('../turnover/src/initialize.R')

ml.merged.files <-  list.files("saved", ".Rdata")


ml.nets <- lapply(file.path("saved", ml.merged.files),
                  function(x){
                      load(x)
                      return(ml.net.all)
                  })

ml.nets.merged <- lapply(file.path("saved", ml.merged.files),
                  function(x){
                      load(x)
                      return(ml.net.merged)
                  })

names(ml.nets) <- names(ml.nets.merged)<-
    gsub(".Rdata", "", ml.merged.files)

page.rank.ml <- lapply(ml.nets,  page_rank)
page.rank.ml.merged <- lapply(ml.nets.merged,  page_rank)
