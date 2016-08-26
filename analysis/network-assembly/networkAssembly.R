rm(list=ls())
library(igraph)
library(bipartite)
library(lme4)
setwd('~/Dropbox/skyIslands/analysis/network-assembly')
source('src/prepNets.R')
source('src/CalcMetrics.R')
source('src/misc.R')
spec <-  read.csv("../data/spec.csv")
N <- 100

## ************************************************************
## create site by year neworks and calculate statistics
## (in terminal)
## ************************************************************
## create pp matrix for each site, year
nets <- break.net(spec)

## ************************************************************
## individuals nulls
## ************************************************************
## nulls <- lapply(nets, vaznull, N=N)
## save(nulls, file='saved/nulls/all.Rdata')

## ************************************************************
## calculate metrics and zscores
## ************************************************************
load(file='saved/nulls/all.Rdata')

mets <- lapply(nets, calc.metric)
null.mets <- rapply(nulls, calc.metric, how="replace")
null.mets <- lapply(null.mets, function(x) do.call(rbind, x))
save(null.mets, file='saved/nullMets.Rdata')

load(file='saved/nullMets.Rdata')

cor.mets <- mapply(function(a, b)
                   cor.metrics(true.stat= a,
                               null.stat= b,
                               N=N),
                   a=mets,
                   b=null.mets,
                   SIMPLIFY=FALSE)

cor.dats <- prep.dat(cor.mets,  spec)
save(cor.dats, file='saved/corMets.Rdata')

## ************************************************
## plotting
load(file='saved/corMets.Rdata')

from.Rockies <- c("JC", "SC", "MM", "PL", "CH")
cor.dats <- cor.dats[match(from.Rockies, cor.dats$Site),]

plot(y=cor.dats$nicheOverlapPol, x= 1:5,
     ylab = "Niche overlap",
     xlab = "Distance from Rockies",
     xaxt="n",
     pch=16, cex=1.5)
axis(1, at=1:5, labels=c("Pecos", "Sandias", "Magdalena", "Pina Leno",
     "Chiricahua"))

plot(y=cor.dats$nicheOverlapPlant, x= 1:5,
     ylab = "Niche overlap",
     xlab = "Distance from Rockies",
     xaxt="n",
     pch=16, cex=1.5)


plot(y=cor.dats$zNODF, x= 1:5,
     ylab = "Nestedness",
     xlab = "Distance from Rockies",
     xaxt="n",
     pch=16, cex=1.5)


plot(y=cor.dats$zmodD, x= 1:5,
     ylab = "Modularity",
     xlab = "Distance from Rockies",
     xaxt="n",
     pch=16, cex=1.5)
