## setwd('~/Dropbox/skyIslands/')
rm(list=ls())
setwd('analysis/network')
source('src/initialize.R')
source('src/vaznull2.R')
net.type <- "YrSR"
## species <- c("Plant", "Pollinator")
species <- c("Pollinator", "Parasite")

source('../turnover/src/initialize.R')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 0){
    N <- as.numeric(args[1])
} else{
    N <- 99
}

## ************************************************************
## calculate metrics and zscores ## beware this takes a while!
## ************************************************************

nets <- nets[apply(sapply(nets, dim) > 3, 2, all))]

mets <- lapply(nets, calcNetworkMetrics,
               N=N)

cor.dats <- prepDat(mets,  spec, net.type=net.type)
save(cor.dats, file=sprintf('saved/corMets_%s.Rdata',
                            paste(species, collapse="")))

## ************************************************************
load(file=sprintf('saved/corMets_%s.Rdata',
                  paste(species, collapse="")))

ys <- c("niche.overlap.LL",
        "niche.overlap.HL",
        "cluster.coefficient.LL",
        "cluster.coefficient.HL",
        "number.of.species.LL",
        "number.of.species.HL",
        "zweighted.NODF",
        "H2")

xvar <- "Lat"

## create formulas for site characteristics
formulas.div <-lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste("scale(DoyPoly1)*scale(Lat)",
                           "scale(DoyPoly2)*scale(Lat)",
                           "Year",
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
     file=file.path(save.path,
                    sprintf('mods/metrics_%s_%s.Rdata',
                            paste(species, collapse=""),
                            xvar)))
