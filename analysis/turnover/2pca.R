## setwd('~/Dropbox/skyIslands')
rm(list=ls())
## This script calculates the pollinator role/network niche
## differences across sites/years
setwd('analysis/turnover')
source('src/calcPca.R')

source('src/initialize.R')


## vector of pca loadings of interest
loadings <- c(1)
metrics <- c("rare.degree",
             "weighted.betweenness",
             "weighted.closeness",
             "niche.overlap",
             "species.strength",
             "d")

sp <- species.roles[species.roles$speciesType == sp.level,]
sp.pca.scores <- calcNetworkPca(species.roles=sp,
                                metrics= metrics,
                                loadings=loadings,
                                nets.by.SR= nets.by.SR)


pcas.species <- split(sp.pca.scores$pcas,
                      sp.pca.scores$pcas$GenusSpecies)

pcas.species <- lapply(pcas.species, calcDiffPcas, geo.dist,
                       nets.by.SR= nets.by.SR)

pcas.species <- pcas.species[!sapply(pcas.species, is.null)]
pcas.species <- do.call(rbind, pcas.species)
rownames(pcas.species) <- NULL


## turnover through space
if(net.type == "YrSR"){
    pcas.diffs.same.year <- pcas.species[apply(pcas.species, 1,
                                        function(x) x["Year1"] ==
                                                    x["Year2"] &
                                                    x["SR1"] ==
                                                    x["SR2"]),]
} else{
    pcas.diffs.same.year <- pcas.species[apply(pcas.species, 1,
                                               function(x){
                                                   x["Year1"] ==
                                                   x["Year2"]}),]
}

## ******************************************************************
## plotting
## ******************************************************************

load(file=sprintf("saved/Beta%s%s.Rdata", net.type,
                  paste(species, collapse="")
                  ))

pcas.beta <- merge(pcas.diffs.same.year, beta.same.year,
                   all.x=TRUE)


xvars <- c("GeoDist", "S", "S_lower.level", "S_higher.level")
xlabs <- c("Geographic distance", "Species turnover",
           paste(species, "turnover"))

modPCATurnover <- function(xvars, beta.same.year){
    this.beta <- pcas.beta
    colnames(this.beta)[colnames(this.beta) == xvars] <- "x"
    forms <- formula(diffPca~scale(x) + (1|Site1) + (1|Site2))
    mod <- do.call(lmer,
                   list(formula=forms,
                        data=this.beta,
                        REML = FALSE))
    eff <- Effect(c("x"), mod)
    return(list(mod=mod, eff=eff))

}
pca.mods <- lapply(xvars, modPCATurnover, beta.same.year)
names(pca.mods) <- xvars

lapply(pca.mods, function(x) summary(x$mod))

source("src/distPlotting.R")
plotPCAs()

save(pcas.diffs.same.year, pcas.species, pca.mods,
     file=sprintf('saved/pcaTurover%s%s%s.Rdata', net.type,
                  paste(species, collapse=""),
                  gsub("[.]", "", sp.level)
                  ))
