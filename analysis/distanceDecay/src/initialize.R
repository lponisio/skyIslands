library(vegan)
library(fields)
source('src/misc.R')
source('src/distDecay.R')

geo <-
  read.csv("../../data/relational/data/relational/tables/geography.csv")

spec <-
  read.csv("../data/spec.csv")

##distance dissimilarity
dist.site <- rdist.earth(cbind(geo$Long, geo$Lat),
                         cbind(geo$Long, geo$Lat))

c.dist <- dist.site[lower.tri(dist.site)]
