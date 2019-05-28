library(vegan)
library(fields)
library(betalink)
source('src/misc.R')

load('../../data/spec.Rdata')
load('../../data/nets.Rdata')

##distance dissimilarity
geo <- unique(spec[, c("Site", "Lat", "Long")])
geo <- geo[!duplicated(geo$Site),]

dist.site <- rdist.earth(cbind(geo$Long, geo$Lat),
                         cbind(geo$Long, geo$Lat))
colnames(dist.site) <- rownames(dist.site) <- geo$Site

c.dist <- dist.site[lower.tri(dist.site)]
