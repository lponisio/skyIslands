## setwd('~/Dropbox/skyIslands')
rm(list=ls())
## This script calculates the pollinator role/network niche
## differences across sites/years
setwd('analysis/turnover')
source('src/calcPca.R')

## net.type <- "YrSR"
net.type <- "Yr"
species <- c("Plant", "Pollinator")
## species <- c("Pollinator", "Parasite")
sp.level <- "lower.level"

source('src/initialize.R')

load(file=sprintf('../../data/splev%s%s.Rdata', net.type,
                  paste(species, collapse="")
                  ))

species.roles <- sp.lev

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
                                               function(x){ x["Year1"] ==
                                                                x["Year2"]}),]
}



save(pcas.diffs.same.year, pcas.species,
     file=sprintf('saved/pcaTurover%s%s%s.Rdata', net.type,
                  paste(species, collapse=""),
                  gsub("[.]", "", sp.level)
                  ))

## ******************************************************************
## plotting
## ******************************************************************

load(file=sprintf("saved/Beta%s%s.Rdata", net.type,
                  paste(species, collapse="")
                  ))

pcas.beta <- merge(pcas.diffs.same.year, beta.same.year,
                   all.x=TRUE)



plotPcaDiff <- function(){
    xvars <- c("GeoDist", "S", "S_lower.level", "S_higher.level")
    xlabs <- c("Geographic distance", "Species turnover",
               paste(species, "turnover"))

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

}

pdf.f(plotPcaDiff,  file=sprintf("figures/pcaDiff%s%s%s.pdf", net.type,
                                 paste(species, collapse=""),
                                 gsub("[.]", "", sp.level)
                                 ),
      height=8, width=8)



