library(vegan)
library(fields)
library(betalink)
source('src/misc.R')
## source('src/distDecay.R')

load('../../data/spec.Rdata')
load('../../data/nets.Rdata')

##distance dissimilarity
geo <- unique(spec[, c("Site", "Lat", "Long")])

dist.site <- rdist.earth(cbind(geo$Long, geo$Lat),
                         cbind(geo$Long, geo$Lat))

c.dist <- dist.site[lower.tri(dist.site)]
