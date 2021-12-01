setwd('C:/Users/rah10/Documents/skyIslands')
## setwd('~/Dropbox (University of Oregon)/skyislands')


setwd("analysis/microbiome")

rm(list=ls())


source("src/init.R")

source("src/misc.R")

ncores <- 1


library(picante)
library(bayesplot)


microbes <- colnames(spec)[grepl("16s:", colnames(spec))]

screened.microbes <- apply(spec, 1, function(x) all(is.na(x[microbes])))

spec.microbes <- spec[!screened.microbes, ]

PD <- apply(spec.microbes[,microbes], 1, function(x){
  this.bee <- x[x > 0]
  this.tree <- prune.sample(t(this.bee), tree.16s)
  pd(t(this.bee), this.tree, include.root = FALSE)
})

PD <- do.call(rbind, PD)

spec.microbes <- cbind(spec.microbes, PD)

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
spec.microbes[, vars] <- apply(spec.microbes[, vars], 2, standardize)


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
fit <- brm(bform, spec.microbes,
           cores=ncores,
           iter = 10^4,
           chains = 2,
           thin=1,
           inits=0,
           control = list(adapt_delta = 0.99))



write.ms.table(fit, "microbes")

save(fit, spec,
     file="saved/microbesFitMod.Rdata")
## dignostic figures
mcmc_trace(fit)
ggsave("figures_diagnostics_microbes.pdf",
       height=11, width=8.5)
Collapse








