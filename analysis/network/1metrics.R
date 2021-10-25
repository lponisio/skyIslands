## setwd('~/Dropbox/skyIslands/')
rm(list=ls())
setwd('analysis/network')
source('src/initialize.R')
source('src/vaznull2.R')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 0){
    N <- as.numeric(args[3])
} else{
    N <- 9
}

## ************************************************************
## calculate metrics and zscores ## beware this takes a while!
## ************************************************************

## nets <- nets[apply(sapply(nets, dim) > 2, 2, all)]

## mets <- lapply(nets, calcNetworkMetrics,
##                N=N)

## cor.dats <- prepDat(mets,  spec, net.type=net.type)

## save(cor.dats, file=sprintf('saved/corMets_%s_%s.Rdata',
##                             paste(species, collapse=""),
##                             net.type))

## ************************************************************
load(file=sprintf('saved/corMets_%s_%s.Rdata',
                            paste(species, collapse=""), net.type))

ys <- c("niche.overlap.LL",
        "niche.overlap.HL",
        "weighted.cluster.coefficient.LL",
        "weighted.cluster.coefficient.HL",
        "mean.number.of.links.LL",
        "mean.number.of.links.HL",
        "number.of.species.LL",
        "number.of.species.HL",
        "zweighted.NODF",
        "zH2")

xvar <- "Lat"

## create formulas for site characteristics
if(species[1] == "Plant"){
    formulas.div <-lapply(ys, function(x) {
        as.formula(paste(x, "~",
                         paste("scale(Doy)*scale(Lat)",
                               "scale(I(Doy^2))*scale(Lat)",
                               ## "scale(Area)",
                               "Year",
                               "(1|Site)",
                               sep="+")))
    })
    mods.div <- lapply(formulas.div, function(x){
        lmer(x, data=cor.dats)
    })
} else if(species[2] == "Parasite"){
    formulas.div <-lapply(ys, function(x) {
        as.formula(paste(x, "~",
                         paste(
                             ## "scale(Area)",
                             "scale(Doy)*scale(Lat)",
                             "scale(I(Doy^2))*scale(Lat)",
                             "(1|Site)",
                             sep="+")))
    })

    mods.div <- lapply(formulas.div, function(x){
        lmer(x, data=cor.dats)
    })

}

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


## are area and latitude correlated?
geo.test <- unique(data.frame(Lat=cor.dats$Lat,
                              Site=cor.dats$Site,
                              Area=cor.dats$Area))

cor.test(geo.test$Lat, log(geo.test$Area))
plot(geo.test$Lat ~log(geo.test$Area))
