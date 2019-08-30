setwd('~/Dropbox/skyIslands/')
rm(list=ls())
setwd('analysis/network')
source('src/initialize.R')
source('src/vaznull2.R')
load('../../data/nets.Rdata')
load('../../data/spec.Rdata')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 0){
    N <- as.numeric(args[1])
} else{
    N <- 99
}

## ************************************************************
## calculate metrics and zscores ## beware this takes a while!
## ************************************************************

## mets <- lapply(nets, calcNetworkMetrics,
##                N=N,   index= c("H2",
##                                "partner diversity",
##                                "functional complementarity",
##                                "links per species",
##                                "weighted NODF",
##                                "NODF",
##                                "number of species",
##                                "niche overlap"))

## cor.dats <- prepDat(mets,  spec)
## save(cor.dats, file='saved/corMets.Rdata')

## ************************************************************
load(file='saved/corMets.Rdata')

ys <- c("functional.complementarity.LL",
        "functional.complementarity.HL",
        "zH2",
        "links.per.species",
        "zweighted.NODF",
        "zmod.met.R",
        "niche.overlap.LL",
        "niche.overlap.HL",
        "partner.diversity.HL",
        "partner.diversity.LL")


## create formulas for site characteristics
formulas.div <-lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste("scale(Lat)",
                           "scale(I(Lat)^2)",
                           "(1|Site)",
                           sep="+")))
})

mods.div <- lapply(formulas.div, function(x){
    lmer(x, data=cor.dats)
})
names(mods.div) <- ys

## results
lapply(mods.div, summary)

## check sig levels with method other than wald CI
lapply(mods.div, anova)

save(mods.div, cor.dats,
     file=file.path(save.path, 'mods/metrics.Rdata'))
