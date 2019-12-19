rm(list=ls())
library(igraph)
library(vegan)
library(fields)
library(fossil)
library(bipartite)
setwd("~/Dropbox/skyIslands/dataPrep")
source('relational/relational_prep.R')
setwd("~/Dropbox/skyIslands/dataPrep")
source('relational/relational_make.R')
setwd("~/Dropbox/skyIslands/dataPrep")
source('relational/make_traditional.R')

setwd("~/Dropbox/skyIslands/dataPrep")
src.dir <- '../../skyIslands_saved/data/relational/relational/traditional/'
spec <-
    read.csv(file.path(src.dir, "specimens-complete.csv"),
             stringsAsFactors=FALSE)

source("src/misc.R")
source("src/prepNets.R")
source("src/specialization.R")

## did not complete full sampling rounds in any of these sites. Was
## just scouting.  can keep UK and SS when more species are IDed
## site.2.drop <- c("JM", "CC", "UK", "SS")
site.2.drop <- c("JM", "CC", "SS", "UK")
spec <- spec[!spec$Site %in% site.2.drop,]
spec <- droplevels(spec)

## get specimen data ready
spec$GenusSpecies <- fix.white.space(paste(spec$Genus,
                          spec$Species,
                          spec$SubSpecies))

spec$PlantGenusSpecies <-  fix.white.space(paste(spec$PlantGenus,
                                          spec$PlantSpecies,
                                          spec$PlantVar,
                                          spec$PlantSubSpecies))

spec$Int <-  fix.white.space(paste(spec$GenusSpecies,
                                   spec$PlantGenusSpecies))
spec$IntGen <-  fix.white.space(paste(spec$Genus,
                                      spec$PlantGenus))

spec$Date <- as.Date(spec$Date, format='%m/%d/%y')
spec$Doy <- as.numeric(strftime(spec$Date, format='%j'))
spec$Year <- as.numeric(format(spec$Date,'%Y'))

## drop non-bees but keep syrphids
spec <- spec[spec$Family %in% c("Andrenidae", "Apidae",
                                "Colletidae", "Halictidae",
                                "Megachilidae", "Syrphidae"),]

## drop non-bees
## spec <- spec[spec$Family %in% c("Andrenidae", "Apidae",
##                                 "Colletidae", "Halictidae",
##                                 "Megachilidae"),]


## drop lasioglossum until we get 2017, 2018 IDs back from Joel
## spec <- spec[spec$Genus != "Lasioglossum",]


## for networks, drop specimens withoutplant IDs
spec <- spec[spec$PlantGenusSpecies != "",]
spec <- spec[spec$GenusSpecies != "",]

## drop the the 2017 sample of PL because it was on fire for other
## sampling rounds and there was basically nothing blooming the first
## round
spec <- spec[!(spec$Site == "PL" & spec$Year == "2017"),]

## calculate orthoganol polynomials for doy
spec$DoyPoly <- poly(spec$Doy, degree=2)
spec$DoyPoly1 <- spec$DoyPoly[,'1']
spec$DoyPoly2 <- spec$DoyPoly[,'2']
spec$DoyPoly <- NULL

## also for Latitude
spec$LatPoly <- poly(spec$Lat, degree=2)
spec$LatPoly1 <- spec$LatPoly[,'1']
spec$LatPoly2 <- spec$LatPoly[,'2']
spec$LatPoly <- NULL

save(spec, file="../data/spec.Rdata")
write.csv(spec, file="../data/spec.csv", row.names=FALSE)

## *******************************************************************
## create a giant plant-pollinator network to calculate specialization
## etc. across all SI
## *******************************************************************
agg.spec <- aggregate(list(abund=spec$GenusSpecies),
                      list(GenusSpecies=spec$GenusSpecies,
                           PlantGenusSpecies=spec$PlantGenusSpecies),
                      length)

nets.all <- samp2site.spp(agg.spec$PlantGenusSpecies,
                          agg.spec$GenusSpecies,
                          agg.spec$abund, FUN=sum)

all.traits <- specieslevel(nets.all)
## calculate rarified plant.pol degree
rare.plants.degree <- apply(nets.all, 1, chao1)
rare.pols.degree <- apply(nets.all, 2, chao1)

traits <- data.frame(GenusSpecies= unlist(sapply(all.traits,
                                                 rownames)),
                     do.call(rbind, all.traits))

traits$r.degree <-  rare.pols.degree[match(traits$GenusSpecies,
                                           names(rare.pols.degree))]
traits$r.degree[is.na(traits$r.degree)] <-
    rare.plants.degree[match(traits$GenusSpecies[is.na(traits$r.degree)],
                             names(rare.plants.degree))]

rownames(traits) <- NULL

write.csv(traits, file='../data/traits.csv')

## *******************************************************************
## create a giant pathogen-pollinator network to calculate
## specialization etc. across all SI
## *******************************************************************
parasites <- c("AspergillusSpp", "AscosphaeraSpp",
                "ApicystisSpp", "CrithidiaExpoeki", "CrithidiaBombi")
                ## "NosemaBombi", "NosemaCeranae")

agg.spec.sub <- spec[!apply(spec[,parasites], 1, function(x) all(is.na(x))),]


## add Nosema when ready
agg.spec.para <- aggregate(agg.spec.sub[, parasites],
                      list(GenusSpecies=agg.spec.sub$GenusSpecies),
                                              sum, na.rm=TRUE)

para.gensp.counts <- table(agg.spec.sub$GenusSpecies)

## proportion of individuals screened
agg.spec.para[, parasites] <- agg.spec.para[, parasites]/
    para.gensp.counts[agg.spec.para$GenusSpecies]

nets.para <- agg.spec.para
rownames(nets.para) <- nets.para$GenusSpecies
nets.para$GenusSpecies <- NULL

all.traits.para <- specieslevel(nets.para)
## calculate rarified plant.pol degree

para.traits <- data.frame(GenusSpecies= unlist(sapply(all.traits.para,
                                                 rownames)),
                     do.call(rbind, all.traits.para))

rownames(para.traits) <- NULL

write.csv(traits, file='../data/parasitetraits.csv')

## *******************************************************************
##  create site, SR, year level networks
## *******************************************************************

## plant-pollinator networks
makeNets(spec, net.type="YrSR")
makeNets(spec, net.type="Yr", mean.by.year=TRUE)



spec.sub <- agg.spec.sub %>%
    select(UniqueID, GenusSpecies, Site, Year, SampleRound, "AspergillusSpp", "AscosphaeraSpp",
                "ApicystisSpp", "CrithidiaExpoeki", "CrithidiaBombi")

prep.para <- spec.sub %>%
    pivot_longer(cols=c("AspergillusSpp", "AscosphaeraSpp",
                        "ApicystisSpp", "CrithidiaExpoeki", "CrithidiaBombi"),
                 names_to = "Parasite", values_to = "count")

prep.para <- as.data.frame(prep.para[prep.para$count == 1,])

makeNets(prep.para, net.type="YrSR", species=c("Pollinator",
                                               "Parasite"),
         lower.level="GenusSpecies",
         higher.level="Parasite")


makeNets(prep.para, net.type="Yr", species=c("Pollinator",
                                               "Parasite"),
         lower.level="GenusSpecies",
         higher.level="Parasite",
         mean.by.year=TRUE)

## *******************************************************************
##  checks
## *******************************************************************


print(paste("Pollinator species", length(unique(spec$GenusSpecies))))
print(paste("Plant species", length(unique(spec$PlantGenusSpecies))))
print(paste("Pollinator genera", length(unique(spec$Genus))))
print(paste("Interactions", length(unique(spec$Int))))
print(paste("Specimens", nrow(spec)))



table(spec$GenusSpecies)

## tab <- table(spec$GenusSpecies, spec$Site)
## tab2 <- table(spec$PlantGenusSpecies, spec$Site)

## table(spec$PlantGenusSpecies)
## table(spec$PlantGenusSpecies, spec$Site)
## table(spec$PlantGenusSpecies, spec$Year)

## table(spec$PlantGenusSpecies, spec$Family)
