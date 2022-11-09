## setwd('C:/Users/na_ma/Dropbox (University of Oregon)/Rotation/skyIslands')
## setwd('~/Dropbox (University of Oregon)/skyislands')
setwd('~/Dropbox (University of Oregon)/skyIslands')

setwd("analysis/microbiome")

rm(list=ls())

source("src/init.R")


source("src/misc.R")

source("../microbiome/src/writeResultsTable.R")

source("../microbiome/src/makeMultiLevelData.R")

ncores <- 1


library(picante)
library(bayesplot)





genus_pd_fit <- function(all_bees, this_genus, num_iter){
  
  microbes <- colnames(all_bees)[grepl("16s:", colnames(all_bees))] 
  
  screened.microbes <- apply(all_bees, 1, function(x) all(is.na(x[microbes])))
  
  spec.microbes <- all_bees[!screened.microbes, ]
  
  genus.microbes <- spec.microbes[spec.microbes$Genus == this_genus, ]
  

PD <- apply(genus.microbes[,microbes], 1, function(x){
  this.bee <- x[x > 0]
  this.tree <- prune.sample(t(this.bee), tree.16s)
  pd(t(this.bee), this.tree, include.root = FALSE)
})

PD <- do.call(rbind, PD)

genus.microbes <- cbind(genus.microbes, PD)

all_bees <- merge(all_bees, genus.microbes, all.x=TRUE)


##copying over code from communityHealthBayes and changing
##parasite for microbiome data

vars <- c("FloralAbundance",
          "FloralDiversity",
          "PollAbundance",
          "PollDiversity",
          "PD",
          "Lat",
          "Elev",
          "Area")

##  center all of the x variables across the datasets
all_bees[, vars] <- apply(all_bees[, vars], 2, standardize)

## will need to modify when we have multiple years
all_bees <- makeDataMultiLevel(all_bees, "Site", "Year")

## **********************************************************
## Model 1.1: formula for forest effects on floral community
## **********************************************************
## flower diversity
formula.flower.div <- formula(FloralDiversity  ~
                                Lat + Area +  (1|Site)
)
## flower abund
formula.flower.abund <- formula(FloralAbundance  ~
                                  Area + (1|Site)
)
## **********************************************************
## Model 1.2: formula for forest effects on bee community
## **********************************************************
## bee diversity
formula.bee.div <- formula(PollDiversity ~
                             FloralAbundance +
                             FloralDiversity +
                             Lat + Area +
                             (1|Site)
)
## bee abund
formula.bee.abund <- formula(PollAbundance ~
                               FloralAbundance +
                               FloralDiversity +
                               Area +
                               (1|Site)
)
## **********************************************************
## Model 1.3: formula for bee community effects on parasitism
## **********************************************************
formula.microbes <- formula(PD ~
                                PollAbundance*FloralDiversity +
                                PollDiversity +
                                FloralAbundance +
                                (1|Site)
)

## **********************************************************
## Community models
## **********************************************************
## convert to brms format
bf.fabund <- bf(formula.flower.abund)
bf.fdiv <- bf(formula.flower.div)
bf.babund <- bf(formula.bee.abund)
bf.bdiv <- bf(formula.bee.div)
bf.microbes <- bf(formula.microbes)

## **********************************************************
## Model 1 community effects on gut microbe phylo distance
## **********************************************************


## full model

bform <- bf.fabund + bf.fdiv + bf.babund + bf.bdiv + bf.microbes +
  set_rescor(FALSE)

## run model
fit <- brm(bform, all_bees,
           cores=ncores,
           iter = num_iter,
           chains = 2,
           thin=1,
           init=0,
           control = list(adapt_delta = 0.99))

return(fit) ## function returns fit
}

#Apis
apis_fit <- genus_pd_fit(spec, 'Apis', 10^3)
write.ms.table(apis_fit, "apis")

#Bombus
bombus_fit <- genus_pd_fit(spec, 'Bombus', 10^3)
write.ms.table(bombus_fit, "bombus")

#Megachile
megachile_fit <- genus_pd_fit(spec, 'Megachile', 10^3)
write.ms.table(megachile_fit, "megachile")

#Anthophora
anthophora_fit <- genus_pd_fit(spec, 'Anthophora', 10^3)
write.ms.table(anthophora_fit, "anthophora")

#write.ms.table(fit, "microbes")


#save(fit, spec,
#     file="saved/microbesFitMod.Rdata")
## dignostic figures
#mcmc_trace(fit)
#ggsave("figures_diagnostics_microbes.pdf",
#       height=11, width=8.5)









