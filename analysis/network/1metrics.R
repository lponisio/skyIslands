rm(list=ls())
setwd("~/")
source("lab_paths.R")

setwd(file.path(local.path, "skyIslands", "analysis", "network"))

source('src/initialize.R')
source('src/vaznull2.R')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 0){
    N <- as.numeric(args[4])
} else{
    N <- 100
}

## ************************************************************
## calculate metrics and zscores ## beware this takes a while!
## ************************************************************

nets <- nets[apply(sapply(nets, dim) >= 2, 2, all)]

mets <- lapply(nets, calcNetworkMetrics,
               N=N)

event <- do.call(rbind, strsplit(names(mets), ".", fixed = TRUE))

cor.dats <- data.frame(
  Site = event[, 1],
  Year = as.integer(event[, 2]),
  SampleRound = event[, 3],
  do.call(rbind, mets),
  row.names = NULL,
  check.names = FALSE
)

save(cor.dats, file = sprintf(
  "saved/corMets_%s_%s_%s_N%d.Rdata",
  paste(species, collapse = ""), net.type, species.group, N
))

## ************************************************************
load(file=sprintf('saved/corMets_%s_%s_%s_N%d.Rdata',
                  paste(species, collapse = ""), net.type, species.group, N
))

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
