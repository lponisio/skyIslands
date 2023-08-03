## setwd("~/Dropbox/skyIslands/")
setwd('~/Dropbox (University of Oregon)/skyIslands/') ## Rebecca wd

rm(list=ls())
setwd("analysis/turnover")
source("src/initialize.R")
source("src/chao.R")
source("src/betaNet.R")
library(ggplot2)
library(lme4)
library(lmerTest)
library(igraph)
library(ggpubr)
library(emmeans)

load("C:/Users/rah10/Dropbox (University of Oregon)/skyIslands/data/networks/microNets.RData")


##adapted from Lauren's 1betalink in skyIslands folder



####################### must use functions from bipartite bc betalink package is depreciated

CH <- spNet_micro$CH
HM <- spNet_micro$HM
JC <- spNet_micro$JC
MM <- spNet_micro$MM
PL <- spNet_micro$PL
SC <- spNet_micro$SC
SM <- spNet_micro$SM

lower.order <- "Microbes"
higher.order <- "Pollinators"


microbe_poll_betalink <- betalinkr_multi(webarray = webs2array(CH, HM, JC, MM, PL, SC, SM),
                                         partitioning="commondenom", binary=TRUE, distofempty='zero', partition.st=TRUE)

#View(microbe_poll_betalink)

colnames(microbe_poll_betalink) <- c("Site1",
                                     "Site2",
                                     "DissimilaritySpeciesComposition",
                                     "OnlySharedLinks",
                                     "WholeNetworkLinks",
                                     "SpeciesTurnoverLinks",
                                     paste("TurnoverAbsence",lower.order,sep=""),
                                     paste("TurnoverAbsence",higher.order,sep=""),
                                     "TurnoverAbsenceBoth"
                                     )






###will need to update LP's function networkBetaDiversity because most of the packages
### are no longer compatible :( 

geo <- unique(spec.net[, c("Site", "Lat", "Long")])
geo <- geo[!duplicated(geo$Site),]

geo.dist <- rdist.earth(cbind(geo$Long, geo$Lat),
                        cbind(geo$Long, geo$Lat))
colnames(geo.dist) <- rownames(geo.dist) <- geo$Site

## add column for geographic distance between sites
microbe_poll_betalink$GeoDist <- apply(microbe_poll_betalink, 1, function(x){
  geo.dist[x["Site1"],  x["Site2"]]
})

#spec.turnover plot this is working!
# 11:38 pm just deciding to copy and paste style it for now lol


# spec.turnover plot

forms <- formula(DissimilaritySpeciesComposition~GeoDist + (1|Site1) + (1|Site2))
#browser()
mod1 <- do.call(lmer,
               list(formula=forms,
                    data=microbe_poll_betalink,
                    REML = FALSE))
summary(mod1)

gr1 <- ref_grid(mod1, cov.keep= c('GeoDist'))
emm1 <- emmeans(gr1, spec= c("GeoDist"), level= 0.95)
emm1

spec.turnover.plot <- ggplot(microbe_poll_betalink, 
                             aes(x=GeoDist, y=DissimilaritySpeciesComposition)) +
  geom_point() + 
  geom_ribbon(data= data.frame(emm1), aes(ymin= lower.CL, ymax= upper.CL, y= NULL), fill= 'grey80', alpha=0.5) +
  geom_line(data= data.frame(emm1), aes(y= emmean)) +
  theme_classic() + 
  labs(x='Geographic Distance (km)', y='Species Turnover') +
  scale_y_continuous(limits=c(0,1))

spec.turnover.plot

# interaction.turnover plot

forms <- formula(WholeNetworkLinks~GeoDist + (1|Site1) + (1|Site2))
#browser()
mod2 <- do.call(lmer,
               list(formula=forms,
                    data=microbe_poll_betalink,
                    REML = FALSE))
summary(mod2)


gr2 <- ref_grid(mod2, cov.keep= c('GeoDist'))
emm2 <- emmeans(gr2, spec= c("GeoDist"), level= 0.95)
emm2


interaction.turnover.plot <- ggplot(microbe_poll_betalink, 
                                    aes(x=GeoDist, y=WholeNetworkLinks)) +
  geom_point() + 
  geom_ribbon(data= data.frame(emm2), aes(ymin= lower.CL, ymax= upper.CL, y= NULL), fill= 'grey80', alpha=0.5) +
  geom_line(data= data.frame(emm2), aes(y= emmean)) +
  theme_classic() + 
  labs(x='Geographic Distance (km)', y='Interaction Turnover') +
  scale_y_continuous(limits=c(0,1))

interaction.turnover.plot

# poll.turnover plot

forms <- formula(TurnoverAbsencePollinators~GeoDist + (1|Site1) + (1|Site2))
#browser()
mod3 <- do.call(lmer,
               list(formula=forms,
                    data=microbe_poll_betalink,
                    REML = FALSE))

summary(mod3)


gr3 <- ref_grid(mod3, cov.keep= c('GeoDist'))
emm3 <- emmeans(gr3, spec= c("GeoDist"), level= 0.95)
emm3



pollinator.turnover.plot <- ggplot(microbe_poll_betalink, 
                                   aes(x=GeoDist, y=TurnoverAbsencePollinators)) +
  geom_point() + 
  geom_ribbon(data= data.frame(emm3), aes(ymin= lower.CL, ymax= upper.CL, y= NULL), fill= 'grey80', alpha=0.5) +
  geom_line(data= data.frame(emm3), aes(y= emmean)) +
  theme_classic() + 
  labs(x='Geographic Distance (km)', y='Species Turnover: Pollinators')+
  scale_y_continuous(limits=c(0,1))

pollinator.turnover.plot

# microbe.turnover plot

forms <- formula(TurnoverAbsenceMicrobes~GeoDist + (1|Site1) + (1|Site2))
#browser()
mod4 <- do.call(lmer,
               list(formula=forms,
                    data=microbe_poll_betalink,
                    REML = FALSE))

summary(mod4)


gr4 <- ref_grid(mod4, cov.keep= c('GeoDist'))
emm4 <- emmeans(gr4, spec= c("GeoDist"), level= 0.95)
emm4


microbe.turnover.plot <- ggplot(microbe_poll_betalink, 
                                aes(x=GeoDist, y=TurnoverAbsenceMicrobes)) +
  geom_point() + 
  geom_ribbon(data= data.frame(emm4), aes(ymin= lower.CL, ymax= upper.CL, y= NULL), fill= 'grey80', alpha=0.5) +
  geom_line(data= data.frame(emm4), aes(y= emmean)) +
  theme_classic() + 
  labs(x='Geographic Distance (km)', y='Species Turnover: Microbes')+
  scale_y_continuous(limits=c(0,1))

microbe.turnover.plot

## int.turnover.sp.comp plot

forms <- formula(SpeciesTurnoverLinks~GeoDist + (1|Site1) + (1|Site2))
#browser()
mod5 <- do.call(lmer,
               list(formula=forms,
                    data=microbe_poll_betalink,
                    REML = FALSE))
summary(mod5)


gr5 <- ref_grid(mod5, cov.keep= c('GeoDist'))
emm5 <- emmeans(gr5, spec= c("GeoDist"), level= 0.95)
emm5


int.turnover.spcomp.plot <- ggplot(microbe_poll_betalink, 
                                   aes(x=GeoDist, y=SpeciesTurnoverLinks)) +
  geom_point() + 
  geom_ribbon(data= data.frame(emm5), aes(ymin= lower.CL, ymax= upper.CL, y= NULL), fill= 'grey80', alpha=0.5) +
  geom_line(data= data.frame(emm5), aes(y= emmean)) +
  theme_classic() + 
  labs(x='Geographic Distance (km)', y='Interaction Turnover: Species Composition') +
  scale_y_continuous(limits=c(0,1))

int.turnover.spcomp.plot

## int.turnover.rewiring plot

forms <- formula(OnlySharedLinks~GeoDist + (1|Site1) + (1|Site2))
#browser()
mod6 <- do.call(lmer,
               list(formula=forms,
                    data=microbe_poll_betalink,
                    REML = FALSE))
summary(mod6)


gr6 <- ref_grid(mod6, cov.keep= c('GeoDist'))
emm6 <- emmeans(gr6, spec= c("GeoDist"), level= 0.95)


int.turnover.rewiring.plot <- ggplot(microbe_poll_betalink, 
                                     aes(x=GeoDist, y=OnlySharedLinks)) +
  geom_point() + 
  geom_ribbon(data= data.frame(emm6), aes(ymin= lower.CL, ymax= upper.CL, y= NULL), fill= 'grey80', alpha=0.5) +
  geom_line(data= data.frame(emm6), aes(y= emmean)) +
  theme_classic() + 
  labs(x='Geographic Distance (km)', y='Interaction Turnover: Rewiring')+
  scale_y_continuous(limits=c(0,1))

int.turnover.rewiring.plot

### facet plots
bee_microbe_turnover <- ggarrange(spec.turnover.plot,
                                  interaction.turnover.plot,
                                  pollinator.turnover.plot,
                                  microbe.turnover.plot,
                                  int.turnover.spcomp.plot,
                                  int.turnover.rewiring.plot,
                                  ncol=2, nrow=3)

bee_microbe_turnover

dir.create("figures")

ggsave(bee_microbe_turnover, file="figures/poll_microbe_betaComponents.pdf",
       height=11, width=8)


##yayyy working up to here! work on function later after ESA

###############################################################################


### having trouble getting function for plotting to work, going to reconstruct plots using ggplot


make_turnover_plot <- function(turnover_type, ylabel, network_betalink){
  
  forms <- formula(turnover_type~GeoDist + (1|Site1) + (1|Site2))
  
  #browser()
  # this error:
  #Error during wrapup: variable lengths differ (found for 'GeoDist')
  mod <- do.call(lmer,
                 list(formula=forms,
                      data=network_betalink,
                      REML = FALSE))
  
  preds <- predict(mod)
  
  
  spec.turnover.plot <- ggplot(network_betalink, 
                               aes(x=GeoDist, y=turnover_type)) +
    geom_point() + 
    geom_smooth(aes(y=preds, x=GeoDist)) +
    theme_classic() + 
    labs(x='Geographic Distance (km)', y=ylabel)
  
  spec.turnover.plot
  
}


yvars <- c("DissimilaritySpeciesComposition",
           "WholeNetworkLinks",
           "TurnoverAbsenceMicrobes",
           "TurnoverAbsencePollinators",
           "SpeciesTurnoverLinks",
           "OnlySharedLinks")
ylabs <- c("Species Turnover", "Interaction Turnover",
           paste("Species Turnover:", lower.order),
           paste("Species Turnover:", higher.order),
           "Interaction Turnover: Species Composition",
           "Interaction Turnover: Rewiring")

make_turnover_plot(turnover_type = yvars[1], ylabel = ylabs[1], network_betalink = microbe_poll_betalink)



for (i in length(yvars)){
  turnover_type <- yvars[i]
  variable_name <- ylabs[i]
  this_plot <- make_turnover_plot(turnover_type = turnover_type,
                                  variable_name = variable_name,
                                  network_betalink = microbe_poll_betalink)
}





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

#yvars <- c("i", "j", "S", "OS", "WN","ST")

ylabs <- c("TurnoverAbsenceMicrobes",
           "TurnoverAbsencePollinators",
           "TurnoverAbsenceBoth",
           "OnlySharedLinks",
           "WholeNetworkLinks",
           "SpeciesTurnoverLinks")

modGeoTurnover <- function(yvars, beta.same.year){
  this.beta <- as.data.frame(beta.same.year)
  y <- yvars
  forms <- formula(y~GeoDist + (1|Site1) + (1|Site2))
  #browser()
  mod <- do.call(lmer,
                 list(formula=forms,
                      data=this.beta,
                      REML = FALSE))
  eff <- Effect(c("GeoDist"), mod)
  return(list(mod=mod, eff=eff))
  
}
geo.mods <- lapply(ylabs, modGeoTurnover, microbe_poll_betalink)
names(geo.mods) <- ylabs

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
