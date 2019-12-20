## setwd("~/Dropbox/skyIslands/")
rm(list=ls())
setwd("analysis/turnover")
net.type <- "Yr"
## net.type <- "Yr"
## species <- c("Plant", "Pollinator")
species <- c("Pollinator", "Parasite")
source("src/initialize.R")
source("src/chao.R")
source("src/betaNet.R")

plants <- unique(spec$PlantGenusSpecies)
pols <- unique(spec$GenusSpecies)
parasites <- c("AspergillusSpp", "AscosphaeraSpp",
               "ApicystisSpp", "CrithidiaExpoeki", "CrithidiaBombi")
## "NosemaBombi", "NosemaCeranae")

if(species[2] == "Pollinator"){
    lower.level  <- plants
    higher.level <- pols
} else if(species[2] == "Parasite"){
    lower.level  <- pols
    higher.level <- parasites
}


beta.net <- networkBetadiversity(nets.graph,
                                  lower.level=lower.level,
                                  higher.level=higher.level,
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

save(beta.same.site, beta.same.year,
     file=sprintf("saved/Beta%s%s.Rdata", net.type,
                      paste(species, collapse="")
                      ))


## ******************************************************************
## plotting
## ******************************************************************

yvars <- c("S", "WN", "S_lower.level", "S_higher.level", "PropST", "OS")


ylabs <- c("Species turnover", "Interaction Turnover",
           paste(species, "turnover"),
           "Interaction turnover due to species turnover",
           "Interaction turnover due to rewiring")

panels <- vector("list", length(yvars))

for(i in 1:length(yvars)){
    this.beta <- beta.same.year
    colnames(this.beta)[colnames(this.beta) == yvars[i]] <- "y"
    panels[[i]] <- ggplot(this.beta,
                          aes(x=GeoDist, y=y,
                              color=Year1)) + geom_point() + geom_smooth(method="lm") +
        labs(x="Geographic Distance", y=ylabs[i]) +
        lims(y=c(0,1))
}

do.call(grid.arrange, panels)


