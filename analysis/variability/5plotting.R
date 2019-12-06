## setwd('~/Dropbox/skyIslands')
rm(list=ls())
setwd('analysis/variability')

load('saved/results/partnerVar_Site.Rdata')
beta.dist.site <- beta.dist
load('saved/results/partnerVar_SiteYear.Rdata')
beta.dist.site.year <- beta.dist
load('saved/results/partnerVar_Year.Rdata')
beta.dist.year <- beta.dist

beta.dist <- rbind(beta.dist.site,
                   beta.dist.site.year,
                   beta.dist.year)

boxplot(beta.dist$dist~beta.dist$Type)

