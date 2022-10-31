rm(list=ls())
library(igraph)
library(vegan)
library(fields)
library(fossil)
library(bipartite)
library(dplyr)
library(RSQLite)
library(tidyr)
library(readr)

dir.bombus <-
    '/Volumes/bombus/Dropbox (University of Oregon)/skyIslands'

## dir.bombus <-
##      '~/Dropbox (University of Oregon)/skyIslands'


setwd(dir.bombus)

source('dataPrep/relational/prep.R')

setwd(dir.bombus)
source('dataPrep/relational/make.R')

setwd(dir.bombus)
source('dataPrep/relational/traditional.R')

## dir.bombus <-
##     '~/Dropbox (University of Oregon)/skyIslands'


dir.bombus <-
    '/Volumes/bombus/Dropbox (University of Oregon)/skyIslands'

setwd(file.path(dir.bombus, "dataPrep"))

src.dir <- '../../skyIslands_saved/data/relational/relational/traditional/'
spec <-
    read.csv(file.path(src.dir, "specimens-complete.csv"),
             stringsAsFactors=FALSE)
bloom <-
    read.csv(file.path(src.dir, "bloom-complete.csv"),
             stringsAsFactors=FALSE)
veg <-
    read.csv(file.path(src.dir, "veg-complete.csv"),
             stringsAsFactors=FALSE)

source("src/misc.R")
source("src/prepNets.R")
source("src/specialization.R")

## did not complete full sampling rounds in any of these sites. Was
## just scouting.  can keep UK and SS when more species are IDed
## site.2.drop <- c("JM", "CC", "SS", "UK")
## spec <- spec[!spec$Site %in% site.2.drop,]
## spec <- droplevels(spec)

## get specimen data ready
spec$GenusSpecies <- fix.white.space(paste(spec$Genus,
                          spec$Species,
                          spec$SubSpecies))

spec$PlantGenusSpecies <-  fix.white.space(paste(spec$PlantGenus,
                                          spec$PlantSpecies,
                                          spec$PlantVar))

spec$Int <-  fix.white.space(paste(spec$GenusSpecies,
                                   spec$PlantGenusSpecies))

spec$Date <- as.Date(spec$Date, format='%Y-%m-%d')
spec$Doy <- as.numeric(strftime(spec$Date, format='%j'))
spec$Year <- as.numeric(format(spec$Date,'%Y'))

## drop non-bees but keep syrphids
spec <- spec[spec$Family %in% c("Andrenidae", "Apidae",
                                 "Colletidae", "Halictidae",
                                 "Megachilidae", "Syrphidae"),]

## for networks, drop specimens withoutplant IDs
spec <- spec[spec$PlantGenusSpecies != "",]
spec <- spec[spec$GenusSpecies != "",]

## drop the the 2017 sample of PL because it was on fire for other
## sampling rounds and there was basically nothing blooming the first
## round, or leave it for phenology?
## spec <- spec[!(spec$Site == "PL" & spec$Year == "2017"),]

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


parasites <- c("AspergillusSpp", "AscosphaeraSpp",
                "ApicystisSpp", "CrithidiaExpoeki", "CrithidiaBombi",
               "NosemaBombi", "NosemaCeranae")

spec[, parasites][is.na(spec[, parasites])] <- 0
spec[, parasites][spec[, parasites] == ""] <- 0
spec[, parasites] <- apply(spec[, parasites], 2,  as.numeric)

spec[spec$Apidae != 1 | is.na(spec$Apidae), parasites] <- NA

spec$ParasiteRichness <- rowSums(spec[, parasites],
                                 na.rm=TRUE)
spec$PossibleParasite <- apply(spec[, parasites], 1,
                               function(x) sum(!is.na(x)))
spec$ParasitePresence <- (spec$ParasiteRichness >= 1)*1

spec[spec$Apidae != 1  | is.na(spec$Apidae), "ParasiteRichness"] <- NA
spec[spec$Apidae != 1  | is.na(spec$Apidae), "ParasitePresence"] <- NA


save(spec, file="../data/spec.Rdata")
write.csv(spec, file="../data/spec.csv", row.names=FALSE)

## ***********************************************************************
## site/species level data
## ***********************************************************************

site.sp <- spec %>%
    group_by(Site, Year, SampleRound, GenusSpecies) %>%
    summarise(Abundance = length(GenusSpecies),
              SpParasitismRate=mean(ParasitePresence, na.rm=TRUE))

site.sum <- spec %>%
    group_by(Site, Year, SampleRound) %>%
    summarise(PollAbundance = length(GenusSpecies),
              PollRichness= length(unique(GenusSpecies)),
              PollDiversity=vegan:::diversity(table(GenusSpecies),
                                              index="shannon"),
              SiteParasitismRate=mean(ParasitePresence, na.rm=TRUE))

hb <- spec[spec$GenusSpecies == "Apis mellifera",]

hb.site.sum <- hb %>%
    group_by(Site, Year, SampleRound) %>%
    summarise(HBAbundance =n(),
              HBSiteParasitismRate=mean(ParasitePresence, na.rm=TRUE))


site.sum  <- merge(site.sum, hb.site.sum, all.x=TRUE)

site.sp.yr <- spec %>%
    group_by(Site, Year, GenusSpecies, Genus) %>%
    summarise(Abundance = length(GenusSpecies))

bombus <- site.sp.yr[site.sp.yr$Genus == "Bombus",]

write.csv(bombus, file='../data/bombus_year_site.csv', row.names=FALSE)
write.csv(site.sp.yr, file='../data/sp_year_site.csv', row.names=FALSE)
write.csv(site.sum, file='../data/sitestats.csv', row.names=FALSE)
write.csv(site.sp, file='../data/spstats.csv', row.names=FALSE)

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

write.csv(traits, file='../data/traits.csv', row.names=FALSE)

## *******************************************************************
## create a giant pathogen-pollinator network to calculate
## specialization etc. across all SI
## *******************************************************************
agg.spec.sub <- spec[spec$Apidae == 1,]

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

write.csv(traits, file='../data/parasitetraits.csv', row.names=FALSE)

## *******************************************************************
##  create site, SR, year level networks
## *******************************************************************

## plant-pollinator networks
makeNets(spec, net.type="YrSR")
makeNets(spec, net.type="Yr", mean.by.year=TRUE)

spec.sub <- agg.spec.sub %>%
    select(UniqueID, GenusSpecies, Site, Year, SampleRound,
           "AspergillusSpp", "AscosphaeraSpp",
           "ApicystisSpp", "CrithidiaExpoeki", "CrithidiaBombi",
            "NosemaBombi", "NosemaCeranae")

prep.para <- spec.sub %>%
    pivot_longer(cols=c("AspergillusSpp", "AscosphaeraSpp",
                        "ApicystisSpp", "CrithidiaExpoeki",
                        "CrithidiaBombi",
                         "NosemaBombi", "NosemaCeranae"),
                 names_to = "Parasite", values_to = "count")
prep.para <- as.data.frame(prep.para)

## prep.para <- as.data.frame(prep.para[prep.para$count == 1,])

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

## table(spec$GenusSpecies)
## tab <- table(spec$GenusSpecies, spec$Site)
## tab2 <- table(spec$PlantGenusSpecies, spec$Site)
## table(spec$PlantGenusSpecies)
## table(spec$PlantGenusSpecies, spec$Site)
## table(spec$PlantGenusSpecies, spec$Year)
## table(spec$PlantGenusSpecies, spec$Family)

## drop non-bees
## spec <- spec[spec$Family %in% c("Andrenidae", "Apidae",
##                                 "Colletidae", "Halictidae",
##                                 "Megachilidae"),]

## spec <- spec[spec$Sex == "f",]

## tab <- table(spec$GenusSpecies, spec$Site, spec$Year)

## all.occ <- apply(tab, 1, sum)

## sp.too.few <- names(all.occ)[all.occ < 5]


## sp.possible <- names(all.occ)[all.occ > 5]

## cleptos <-  c("Sphecodes sp. a", "Sphecodes sp. b", "Triepeolus sp. a")

## sp.possible <- sp.possible[!sp.possible %in% cleptos]

## sp.sub <- spec[spec$GenusSpecies %in% sp.possible,]

## tab.sub <- sp.sub  %>%
##     group_by(GenusSpecies, Site, Year, SampleRound) %>%
##     summarise(n = n())


##  tab.sub$n30p <- ceiling(tab.sub$n*0.3)
##  tab.sub$n20p <- ceiling(tab.sub$n*0.2)

## sum(tab.sub$n30p)
## sum(tab.sub$n20p)

## tab.sub.cutoff <- tab.sub

## tab.sub.cutoff$n[tab.sub.cutoff$n > 5] <- 5
## sum(tab.sub.cutoff$n)

## *******************************************************************
##  Veg and bloom
## *******************************************************************

veg$PlantGenusSpecies <-  fix.white.space(paste(veg$PlantGenus,
                                          veg$PlantSpecies,
                                          veg$PlantVar))
veg <- veg[veg$PlantGenusSpecies != "",]

bloom$PlantGenusSpecies <-  fix.white.space(paste(bloom$PlantGenus,
                                          bloom$PlantSpecies,
                                          bloom$PlantVar))
bloom <- bloom[bloom$PlantGenusSpecies != "",]

sort(unique(veg$PlantGenusSpecies))
sort(unique(bloom$PlantGenusSpecies))

veg$Date <- as.Date(veg$Date,  "%m/%d/%y")
veg$Year <- format(veg$Date,  "%Y")

bloom$Date <- as.Date(bloom$Date,  "%m/%d/%y")
bloom$Year <- format(bloom$Date,  "%Y")

## blooming flowers
veg.blooming <- veg[veg$BloomStatus != "not blooming" &
                    veg$BloomStatus != "" &
                    veg$NumBlooms != "" &
                    veg$NumBlooms != "0",]

veg.blooming$NumBloomsNum <- 0
veg.blooming$NumBloomsNum[veg.blooming$NumBlooms == "<10"] <- 10
veg.blooming$NumBloomsNum[veg.blooming$NumBlooms == "10-100"] <- 100
veg.blooming$NumBloomsNum[veg.blooming$NumBlooms == "100-1000"] <- 1000

veg.bloom.sum.sp <- veg.blooming %>%
    group_by(Site, Year, PlantGenusSpecies, SampleRound) %>%
    summarise(FloralAbundance = sum(NumBloomsNum))

veg.bloom.sum.sp$SiteDate <- paste(veg.bloom.sum.sp$Site,
                               veg.bloom.sum.sp$SampleRound)

floral.div <- tapply(veg.bloom.sum.sp$FloralAbundance,
                     veg.bloom.sum.sp$SampleRound,
       vegan:::diversity)

veg.bloom.sum.sp <- veg.blooming %>%
    group_by(Site, Year, SampleRound, PlantGenusSpecies) %>%
    summarise(SpFloralAbundance = sum(NumBloomsNum))

veg.bloom.sum <- veg.bloom.sum.sp %>%
    group_by(Site, Year, SampleRound) %>%
    summarise(FloralAbundance = sum(SpFloralAbundance),
              FloralRichness= length(unique(PlantGenusSpecies)),
              FloralDiversity=vegan:::diversity(SpFloralAbundance,
                                                 index="shannon"))



veg.bloom.sum$SiteSR <- paste(veg.bloom.sum$Site,
                               veg.bloom.sum$SampleRound)

veg.bloom.sum$FloralDiv <- floral.div[match(veg.bloom.sum$SiteSR,
                                            names(floral.div))]

veg.bloom.sum$SiteSR <- NULL

write.csv(veg.bloom.sum, file="../data/veg.csv", row.names=FALSE)


## *******************************************************************
## checking data between years
## *******************************************************************

to.sample <- c("CH", "HM", "JC", "MM", "PL", "RP", "SC", "SM")

veg.prep <- veg[veg$Site %in% to.sample,]
veg.split <- split(veg.prep, veg.prep$Site)

spec.veg <- spec[spec$Site %in% to.sample,]
spec.split <- split(spec.veg, spec.veg$Site)

checkveg <- function(x){
    print(unique(x$Site))
    years <- tapply(x$PlantGenusSpecies, x$Year, unique)
    if(length(years) == 2){
        badmatch12 <- years[[1]][!years[[1]] %in% years[[2]]]
        badmatch21 <- years[[2]][!years[[2]] %in% years[[1]]]
        return(list(badmatch12=badmatch12,
                    badmatch21=badmatch21))
    } else if(length(years) == 3){
        badmatch12 <- years[[1]][!years[[1]] %in% years[[2]]]
        badmatch21 <- years[[2]][!years[[2]] %in% years[[1]]]

        badmatch13 <- years[[1]][!years[[1]] %in% years[[3]]]
        badmatch31 <- years[[3]][!years[[3]] %in% years[[1]]]

        badmatch23 <- years[[2]][!years[[2]] %in% years[[3]]]
        badmatch32 <- years[[3]][!years[[3]] %in% years[[2]]]

        return(list(badmatch12=badmatch12,
                    badmatch21=badmatch21,
                    badmatch13= badmatch13,
                    badmatch31=badmatch31,
                    badmatch23=badmatch23,
                    badmatch32=badmatch32))
    } else{
        return(NULL)
    }
}

badmatches.veg <- lapply(veg.split, checkveg)

badmatches.spec <- lapply(spec.split, checkveg)


all.plant.sp <- rbind(data.frame(Site=veg.prep$Site,
                                 PlantGenusSpecies=
                                     veg.prep$PlantGenusSpecies),
                      data.frame(Site=spec.veg$Site,
                                 PlantGenusSpecies=
                                     spec.veg$PlantGenusSpecies))

all.plant.sp <- unique(all.plant.sp)
all.plant.sp <- sort(all.plant.sp)

write.csv(all.plant.sp, row.names=FALSE,
          file="../data/vegbySite.csv")
