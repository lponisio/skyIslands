## setwd("~/Dropbox/skyIslands/")
rm(list=ls())
setwd("analysis/turnover")
## net.type <- "YrSR"
net.type <- "Yr"
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
if(net.type == "YrSR"){
    nets.by.SR  <- TRUE
} else {
    nets.by.SR  <- FALSE
}

beta.net <- networkBetadiversity(nets.graph,
                                 lower.level=lower.level,
                                 higher.level=higher.level,
                                 geo.dist=geo.dist,
                                 nets.by.SR=nets.by.SR)

## Bs: Dissimilarity in the species composition of communities
## Bwn: Dissimilarity of interactions
## Bos: Dissimilarity of interactions established between species
## common to both realisations
## Bst: Dissimilarity of interactions due to species turnover
## ProbST: Bst/wn: Dissimilarity of interactions due to species turnover

## turnover though time
beta.same.site <- beta.net[apply(beta.net, 1,
                                 function(x) x["Site1"] ==
                                             x["Site2"] &
                                             x["Year1"] !=
                                             x["Year2"]),]

## turnover through space
beta.same.year <- beta.net[apply(beta.net, 1,
                                 function(x) x["Year1"] ==
                                             x["Year2"] &
                                             x["Site1"] !=
                                             x["Site2"]),]

## ## turnover though time within a year
beta.same.site.year <- beta.net[apply(beta.net, 1,
                                      function(x) x["Site1"] ==
                                                  x["Site2"] &
                                                  x["Year1"] ==
                                                  x["Year2"]),]



save(beta.same.site, beta.same.year, same.site.year,
     file=sprintf("saved/Beta%s%s.Rdata", net.type,
                  paste(species, collapse="")
                  ))


## ******************************************************************
## plotting
## ******************************************************************

yvars <- c("S", "WN", "S_lower.level", "S_higher.level", "PropST",
           "OS")

ylabs <- c("Species Turnover", "Interaction Turnover",
           paste("Species Turnover:", species),
           "Interaction Turnover: Species Composition",
           "Interaction Turnover: Rewiring")

panels <- vector("list", length(yvars))

plotGeoTurnover <- function(){
    for(i in 1:length(yvars)){
        this.beta <- beta.same.year
        colnames(this.beta)[colnames(this.beta) == yvars[i]] <- "y"
        panels[[i]] <- ggplot(this.beta,
                              aes(x=GeoDist, y=y)) +
            geom_point() + geom_smooth(method="lm") +
            labs(x="Geographic Distance", y=ylabs[i]) +
            lims(y=c(0,1))
    }

    do.call(grid.arrange, panels)
}

pdf.f(plotGeoTurnover,  file=sprintf("figures/GeoBeta%s%s.pdf", net.type,
                                     paste(species, collapse="")
                                     ),
      height=9, width=8)

plotYrTurnover <- function(){
    for(i in 1:length(yvars)){
        this.beta <- beta.same.site
        colnames(this.beta)[colnames(this.beta) == yvars[i]] <- "y"

        panels[[i]] <- ggplot(this.beta) +
            geom_boxplot(aes(factor(YrDist), y=y)) +
            scale_x_discrete(labels = c("1", "5", "6")) +
            labs(x="Years between surveys", y=ylabs[i]) +
            lims(y=c(0,1))
    }

    do.call(grid.arrange, panels)
}

pdf.f(plotYrTurnover,  file=sprintf("figures/YrBeta%s%s.pdf", net.type,
                                    paste(species, collapse="")
                                    ),
      height=9, width=8)

if(nets.by.SR){
    plotSRTurnover <- function(){
        for(i in 1:length(yvars)){
            this.beta <- beta.same.site.year
            colnames(this.beta)[colnames(this.beta) == yvars[i]] <- "y"
            panels[[i]] <- ggplot(this.beta,
                                  aes(x=SRDist, y=y)) +
                geom_point() + geom_smooth(method="lm") +
                labs(x="Days", y=ylabs[i]) +
                lims(y=c(0,1))
        }

        do.call(grid.arrange, panels)
    }

    pdf.f(plotSRTurnover,  file=sprintf("figures/SRBeta%s%s.pdf", net.type,
                                        paste(species, collapse="")
                                        ),
          height=9, width=8)
}
