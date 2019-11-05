## setwd('~/Dropbox/skyIslands/')
rm(list=ls())
setwd('analysis/network')
source('src/initialize.R')
load('../../data/nets.Rdata')
load('../../data/spec.Rdata')

## either "abund" or "degree"
extinction.method <- "abund"

## **********************************************************
## robustness
## **********************************************************
## simulation plant extinction

res <- simExtinction(nets, extinction.method, spec)

save(res, file=file.path(save.path,
            sprintf('resilience_%s.Rdata', extinction.method)))

mod.latitude <- lmer(Robustness ~ Lat
             + (1|Site) + (1|Year),
             data=res)
summary(mod.latitude)
save(mod.latitude, file=file.path(save.path,
            sprintf('mods/resilience_status_%s.Rdata', extinction.method)))

