rm(list=ls())
setwd('~/Dropbox/skyIslands/analysis/networkLevel')
source('src/initialize.R')
load('../../data/nets.Rdata')
load('../../data/spec.Rdata')
N <- 99

## ## ************************************************************
## ## calculate metrics and zscores
## ## ************************************************************
mets <- lapply(nets, network.metrics, N)

cor.dats <- prep.dat(mets,  spec)

cor.dats$tot.rich <- cor.dats$number.of.species.LL +
  cor.dats$number.of.species.HL

## ************************************************************
## degree dists
## ************************************************************
dd.HL <- lapply(nets, colSums)
dd.LL <- lapply(nets, rowSums)

dats.degree <- data.frame(degree=do.call(c, c(dd.HL, dd.LL)),
                          site.type = c(rep(names(dd.HL),
                                            sapply(dd.HL, length)),
                                        rep(names(dd.LL),
                                            sapply(dd.LL, length))),
                          species= rep(c("pollinator", "plant"),
                                       c(length(unlist(dd.HL)),
                                         length(unlist(dd.LL)))))

dats.degree$GenusSpecies <- sapply(strsplit(rownames(dats.degree), "[.]"),
                                   function(x) x[3])

dats.degree$Site <-
    sapply(strsplit(as.character(dats.degree$site.type),
                    "[.]"), function(x) x[1])
dats.degree$Year <-
    sapply(strsplit(as.character(dats.degree$site.type),
                    "[.]"), function(x) x[2])
rownames(dats.degree) <- NULL

## mean degree for each site

mean.degree <- aggregate(list(mean.degree=dats.degree$degree),
                         list(Site=dats.degree$Site,
                              Year=dats.degree$Year,
                              SpeciesType=dats.degree$species),
                         function(x) c(sd(x), mean(x)))
mean.degree$mean <- mean.degree$mean.degree[,1]
mean.degree$sd <- mean.degree$mean.degree[,2]
mean.degree$mean.degree <- NULL

mean.degree.pol <- mean.degree[mean.degree$SpeciesType == "pollinator",]
mean.degree.plant <- mean.degree[mean.degree$SpeciesType == "plant",]

cor.dats <- merge(cor.dats, mean.degree.pol)
colnames(cor.dats)[colnames(cor.dats) %in% c("mean", "sd")] <-
    c("pol.mean.degree", "pol.sd.degree")
cor.dats$SpeciesType <- NULL

cor.dats <- merge(cor.dats, mean.degree.plant)

colnames(cor.dats)[colnames(cor.dats) %in% c("mean", "sd")] <-
    c("plant.mean.degree", "pland.sd.degree")

geo <- unique(spec[, c("Site", "Lat", "Long")])
geo <- geo[!duplicated(geo$Site),]

cor.dats <- merge(cor.dats, geo)

save(cor.dats, file='saved/corMets.Rdata')

