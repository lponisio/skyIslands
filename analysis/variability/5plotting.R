## setwd('~/Dropbox/skyIslands')
rm(list=ls())
setwd('analysis/variability')
source("src/misc.R")

load('saved/results/partnerVar_Site.Rdata')
beta.dist.site <- beta.dist
load('saved/results/partnerVar_SiteYear.Rdata')
beta.dist.site.year <- beta.dist
load('saved/results/partnerVar_Year.Rdata')
beta.dist.year <- beta.dist

beta.dist <- rbind(beta.dist.site,
                   beta.dist.site.year,
                   beta.dist.year)

f <- function(){
    boxplot(beta.dist$dist~beta.dist$Type,
            names=c("Btwn sites", "Btwn surveys", "Btwn years"),
            xlab="", ylab="Interaction beta diversity")
}

pdf.f(f, file= "figures/IntBeta.pdf")
