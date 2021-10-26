rm(list=ls())
library(vegan)
library(fields)
library(igraph)
library(linkcomm)
library(bipartite)
setwd("~/Dropbox/Sky Islands/Analysis/")
source('functions/samp2site_spp.R')
source('~/Dropbox/hedgerow/network/assembly/src/calc_metrics.R')

geo <-
  read.csv("../Data/Relational/data/relational/tables/geography.csv")

spec <-
  read.csv("data/spec.csv")

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}


## plants as sites, polinators as species
prep.comm <- aggregate(spec$GenSp,
                       list(site=spec$PlantGenSp,
                            sp=spec$GenSp), length)

comm <- samp2site.spp(prep.comm$site, prep.comm$sp, prep.comm$x)
comm.dis <- as.matrix(vegdist(comm, "chao", diag= TRUE))

## mds2 <- wcmdscale(comm.dis, k=2)
## plot(mds2[,1], mds2[,2], type = "n", xlab = "", ylab = "",
##      axes = FALSE, main = "wcmdscale (vegan)")
## text(mds2[,1], mds2[,2], rownames(comm.dis), cex = 0.5)

##calculate specialization
d <- specieslevel(comm, index="d")

## pollinators as sites, plants as species
prep.comm <- aggregate(spec$GenSp,
                       list(site=spec$GenSp,
                            sp=spec$PlantGenSp), length)

comm <- samp2site.spp(prep.comm$site, prep.comm$sp, prep.comm$x)
comm.dis <- as.matrix(vegdist(comm, "chao", diag= TRUE))

## pollinators as sites, sites as species

prep.comm <- aggregate(spec$GenSp,
                       list(site=spec$GenSp,
                            sp=spec$Site), length)

comm <- samp2site.spp(prep.comm$site, prep.comm$sp, prep.comm$x)
comm.dis <- as.matrix(vegdist(comm, "chao", diag= TRUE))

evals <- eigen(comm.dis, only.values=TRUE, symmetric=TRUE)

evals <- matrix(evals$values, ncol=1)
rownames(evals) <- rownames(comm.dis)

lm.evals <- lm(evals ~ d$'higher level'$d)

summary(lm.evals)

plot(evals ~ d$'higher level'$d)


## sites as sites, pollinators as species

prep.comm <- aggregate(spec$GenSp,
                       list(site=spec$Site,
                            sp=spec$GenSp), length)

comm <- samp2site.spp(prep.comm$site, prep.comm$sp, prep.comm$x)
comm.dis <- as.matrix(vegdist(comm, "chao", diag= TRUE))

mds1 <- metaMDS(comm)

mds2 <- decorana(comm)





