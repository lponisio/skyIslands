setwd('/Volumes/bombus/Dropbox (University of Oregon)/skyislands')

## setwd('~/Dropbox (University of Oregon)/skyislands')

setwd("analysis/microbiome")
rm(list=ls())

source("src/init.R")
source("src/misc.R")

library(geiger)
library(picante)
library(adiv)

microbes <- colnames(spec)[grepl("16s:", colnames(spec))]

screened.microbes <- apply(spec, 1, function(x) all(is.na(x[microbes])))
spec.microbes <- spec[!screened.microbes, ]

PD <- apply(spec.microbes[,microbes], 1, function(x){
    this.bee <- x[x > 0]
    this.tree <- prune.sample(t(this.bee), tree.16s)
    pd(t(this.bee), this.tree, include.root = FALSE)
})

PD <- do.call(rbind, PD)

spec.microbes <- cbind(spec.microbes, PD)




all.phylo.dist <- apply(spec.microbes[,microbes], 1, function(x){
    this.bee <- x[x > 0]*1000
    browser()
    this.tree <- prune.sample(t(this.bee), tree.16s)
    out <- this.evo.div(this.tree, t(this.bee), tol = 1e-20)
})

