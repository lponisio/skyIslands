## setwd('~/Dropbox/skyislands')
rm(list=ls())
setwd('analysis/variability')
source('src/initialize.R')

net.types <- c("Year", "Site")

## when net.type == "Year", beta diversity is measured for each
## species to examine turnover between sites within a year

## when net.type == "Site", beta diversity is measured for each
## species to examine turnover between years within a site

betas <- list()
for(net in net.types){
    load(file=sprintf("saved/results/partnerVar_%s.Rdata",
                      net))
    beta.dist$BetaType <- net
    betas[[net]] <- beta.dist
}


betas <- do.call(rbind, betas)
rownames(betas) <- NULL

mod.beta <- lmer(dist~BetaType + (1|Site) +
         (1|Year) + (1|GenusSpecies),
     data=betas)
