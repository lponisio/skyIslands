## setwd("~/Dropbox/skyIslands/")
rm(list=ls())
setwd("analysis/turnover")
source("src/initialize.R")
source("src/chao.R")
source("src/betaNet.R")

plants <- unique(spec$PlantGenusSpecies)
pols <- unique(spec$GenusSpecies)

nets.graph <- prepare_networks(nets)

beta.net <- network_betadiversity(nets.graph,
                                  plants=plants, pols=pols,
                                  geo.dist=geo.dist)

## Bs: Dissimilarity in the species composition of communities
## Bwn: Dissimilarity of interactions
## Bos: Dissimilarity of interactions established between species
## common to both realisations
## Bst: Dissimilarity of interactions due to species turnover
## ProbST: Bst/wn: Dissimilarity of interactions due to species turnover

## turnover though time
beta.same.site <- beta.net[apply(beta.net, 1,
                                 function(x) x["Site1"] ==
                                             x["Site2"]),]

## turnover through space
beta.same.year <- beta.net[apply(beta.net, 1,
                                 function(x) x["Year1"] ==
                                             x["Year2"] &
                                             x["Site1"] !=
                                             x["Site2"]),]

save(beta.same.site, beta.same.year, file="saved/SpIntTurnover.Rdata")


yvars <- c("S", "WN", "S_Plants", "S_Pols", "PropST", "OS")

ylabs <- c("Species turnover", "Interaction Turnover",
           "Plant turnover",
           "Pollinator turnover",
           "Interaction turnover due to species turnover",
           "Interaction turnover due to rewiring")

panels <- vector("list", length(yvars))

for(i in 1:length(yvars)){
    this.beta <- beta.same.year
    colnames(this.beta)[colnames(this.beta) == yvars[i]] <- "y"
    panels[[i]] <- ggplot(this.beta,
       aes(x=GeoDist, y=y,
           color=Year1)) + geom_point() + geom_smooth() +
    labs(x="Geographic Distance", y=ylabs[i]) +
    lims(y=c(0,1))
}

do.call(grid.arrange, panels)


