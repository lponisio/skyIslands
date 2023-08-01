## setwd("~/Dropbox/skyIslands/")
setwd('~/Dropbox (University of Oregon)/skyIslands/') ## Rebecca wd
#rm(list=ls())
setwd("analysis/turnover")
source("src/initialize.R")
source("src/chao.R")
source("src/betaNet.R")

library(lme4)
library(igraph)

#need to src microNets.R


##adapted from Lauren's 1betalink in skyIslands folder



####################### must use functions from bipartite bc betalink package is depreciated

CH <- spNet_micro$CH
HM <- spNet_micro$HM
JC <- spNet_micro$JC
MM <- spNet_micro$MM
PL <- spNet_micro$PL
SC <- spNet_micro$SC
SM <- spNet_micro$SM


microbe_poll_betalink <- betalinkr_multi(webarray = webs2array(CH, HM, JC, MM, PL, SC, SM),
                                         partitioning="poisot", binary=TRUE, distofempty='na')

microbe_poll_betalink

###will need to update LP's function networkBetaDiversity because most of the packages
### are no longer compatible :( 

geo <- unique(spec.net[, c("Site", "Lat", "Long")])
geo <- geo[!duplicated(geo$Site),]

geo.dist <- rdist.earth(cbind(geo$Long, geo$Lat),
                        cbind(geo$Long, geo$Lat))
colnames(geo.dist) <- rownames(geo.dist) <- geo$Site

## add column for geographic distance between sites
microbe_poll_betalink$GeoDist <- apply(microbe_poll_betalink, 1, function(x){
  geo.dist[x["i"],  x["j"]]
})

# i = Site1
# j = Site2
# S = dissimilarity in species composition
# OS = dissimilarity explained by rewiring among shared species (only shared)
# WN = dissimilarity between two networks (whole network)
# ST = dissimilarity explained by difference in species community composition (species turnover links)


# for identical results to poisot 2012 betalink use the following settings:
#
#partitioning="poisot", function.dist="betadiver", distofempty="na" and binary=TRUE
# including the function.dist induces a weird error.... need to figure out still if we want to use this method


## ******************************************************************
## calculate different breakdowns of turnover
## ******************************************************************


# beta.net <- networkBetadiversity(microbe_igraph_list,
#                                  lower.level=pols,
#                                  higher.level=microbes,
#                                  geo.dist=geo.dist,
#                                  nets.by.SR=nets.by.SR)
# 
# ## beta.net is formed BUT -- OS, ST, S.lower level, S. higher level, prop ST are induced NaNs
# # need to examine function more to understand why this is happening
# # ALSO check geo.dist to make sure the distance between the same sites are zero,
# # it looks like rn JC, SM, HM, and MM are not 0 distance between themselves :( WHYYYYY
# 
# 
# ## turnover though time
# # beta.same.site <- beta.net[apply(beta.net, 1,
# #                                  function(x) x["Site1"] ==
# #                                              x["Site2"] &
# #                                              x["Year1"] !=
# #                                              x["Year2"]),]
# ## turnover through space
# beta.same.year <- beta.net[apply(beta.net, 1,
#                                  function(x) x["Year1"] ==
#                                              x["Year2"] &
#                                              x["Site1"] !=
#                                              x["Site2"]),]
# ## ## turnover though time within a year
# beta.same.site.year <- beta.net[apply(beta.net, 1,
#                                       function(x) x["Site1"] ==
#                                                   x["Site2"] &
#                                                   x["Year1"] ==
#                                                   x["Year2"]),]

#### still only one year of microbe data so not an issue yet

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
