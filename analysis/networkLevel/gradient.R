rm(list=ls())
setwd('~/Dropbox/skyIslands/analysis/networkLevel')
source('src/initialize.R')
load('../data/networks/nets.Rdata')
load('../data/spec/spec.Rdata')
N <- 99

## ## ************************************************************
## ## calculate metrics and zscores
## ## ************************************************************
mets <- lapply(nets, network.metrics, N)

cor.dats <- prep.dat(mets,  spec)

cor.dats$tot.rich <- cor.dats$number.of.species.LL +
  cor.dats$number.of.species.HL

save(cor.dats, file='saved/corMets.Rdata')

## ************************************************************
## effect of ??
## ************************************************************
load(file='saved/corMets.Rdata')

