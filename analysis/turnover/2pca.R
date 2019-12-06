## This script calculates the pollinator role/network niche
## differences across sites/years
## setwd('~/Dropbox/skyIslands')
rm(list=ls())
setwd('analysis/turnover')
source('src/initialize.R')
source('src/calcPca.R')
source('src/calcSpec.R')


species.roles <- calcSpec(nets, spec, dist.metric="chao")

## vector of pca loadings of interest
loadings <- c(1)
metrics <- c("rare.degree",
             "weighted.betweenness",
             "weighted.closeness",
             "niche.overlap",
             "species.strength",
             "d")

pol <- species.roles[species.roles$speciesType == "pollinator",]
pol.pca.scores <- calcNetworkPca (species.roles=pol,
                                  metrics= metrics,
                                  loadings=loadings)


pcas.species <- split(pol.pca.scores$pcas,
                      pol.pca.scores$pcas$GenusSpecies)

pcas.species <- lapply(pcas.species, calcDiffPcas, geo.dist)
pcas.species <- pcas.species[!sapply(pcas.species, is.null)]
pcas.species <- do.call(rbind, pcas.species)
rownames(pcas.species) <- NULL


## turnover through space
pcas.diffs.same.year <- pcas.species[apply(pcas.species, 1,
                                           function(x) x["Year1"] ==
                                                       x["Year2"] &
                                                       x["SR1"] ==
                                                       x["SR2"]),]

save(pcas.diffs.same.year, file="saved/PcaTurnover.Rdata")


load(file="saved/SpIntTurnover.Rdata")

pcas.beta <- merge(pcas.diffs.same.year, beta.same.year,
                   all.x=TRUE)

## ******************************************************************
## plotting
## ******************************************************************

xvars <- c("GeoDist", "S", "S_Plants", "S_Pols")
xlabs <- c("Geographic distance", "Species turnover",
           "Plant turnover",
           "Pollinator turnover")

panels <- vector("list", length(xvars))

for(i in 1:length(xvars)){
    this.beta <- pcas.beta
    colnames(this.beta)[colnames(this.beta) == xvars[i]] <- "x"
    panels[[i]] <- ggplot(this.beta,
       aes(x=x, y=diffPca,
           color=GenusSpecies, pch=Year1)) + geom_point()  +
    labs(y="Change in network position", x=xlabs[i])  +
        theme(legend.position = "none")
}

do.call(grid.arrange, panels)


