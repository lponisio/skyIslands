library(vegan)
library(fields)
library(betalink)
library(ggplot2)
library(gridExtra)

source('src/misc.R')

load('../../data/spec.Rdata')
load('../../data/netsYrSR.Rdata')


##distance dissimilarity
geo <- unique(spec[, c("Site", "Lat", "Long")])
geo <- geo[!duplicated(geo$Site),]

geo.dist <- rdist.earth(cbind(geo$Long, geo$Lat),
                         cbind(geo$Long, geo$Lat))
colnames(geo.dist) <- rownames(geo.dist) <- geo$Site
