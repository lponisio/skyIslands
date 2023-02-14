## setwd("~/Dropbox/skyIslands/")
setwd('~/Dropbox (University of Oregon)/skyIslands/') ## Rebecca wd
rm(list=ls())
setwd("analysis/turnover")
source("src/initialize.R")
source("src/chao.R")
source("src/betaNet.R")



##adapted from Lauren's 1betalink in skyIslands folder

## need to make cooccurrence network of bees and microbes at each site
## then use betalinkPP function to find the betadiversity between each site
## betalinkPP <- function (n1, n2, bf = B01, lower.level, higher.level)
## do we need to do this for each combination?
## what is lower.level and higher.level?
##      i think higher level will be bees and lower level will be microbes?
## what is nets.graph?

## ******************************************************************
## calculate different breakdowns of turnover
## ******************************************************************


beta.net <- networkBetadiversity(nets.graph,
                                 lower.level=lower.level,
                                 higher.level=higher.level,
                                 geo.dist=geo.dist,
                                 nets.by.SR=nets.by.SR)
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

## ******************************************************************
## plotting/models
## ******************************************************************
## Bs: Dissimilarity in the species composition of communities
## Bwn: Dissimilarity of interactions
## Bos: Dissimilarity of interactions established between species
## common to both realisations
## Bst: Dissimilarity of interactions due to species turnover
## ProbST: Bst/wn: Dissimilarity of interactions due to species turnover

yvars <- c("S", "WN", "S_lower.level", "S_higher.level", "PropST",
           "OS")

ylabs <- c("Species Turnover", "Interaction Turnover",
           paste("Species Turnover:", species),
           "Interaction Turnover: Species Composition",
           "Interaction Turnover: Rewiring")


modGeoTurnover <- function(yvars, beta.same.year){
    this.beta <- beta.same.year
    colnames(this.beta)[colnames(this.beta) == yvars] <- "y"
    forms <- formula(y~scale(GeoDist) + (1|Site1) + (1|Site2))
    mod <- do.call(lmer,
                   list(formula=forms,
                        data=this.beta,
                        REML = FALSE))
    eff <- Effect(c("GeoDist"), mod)
    return(list(mod=mod, eff=eff))

}
geo.mods <- lapply(yvars, modGeoTurnover, beta.same.year)
names(geo.mods) <- yvars

lapply(geo.mods, function(x) summary(x$mod))

source("src/distPlotting.R")
plotDists()


save(beta.same.site, beta.same.year, beta.same.site.year, geo.mods,
     file=sprintf("saved/Beta%s%s.Rdata", net.type,
                  paste(species, collapse="")
                  ))
