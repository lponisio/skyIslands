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

## dir.bombus <-
##     '/Volumes/bombus/Dropbox (University of Oregon)/skyIslands'
dir.bombus <-
     '~/Dropbox (University of Oregon)/skyIslands'

## *****************************************************************
## create relational database, add species IDs
## *****************************************************************

setwd(dir.bombus)
source('dataPrep/relational/prep.R')

setwd(dir.bombus)
source('dataPrep/relational/make.R')

setwd(dir.bombus)
source('dataPrep/relational/traditional.R')

## *****************************************************************
## prep specimen data
## *****************************************************************
dir.bombus <-
    '~/Dropbox (University of Oregon)/skyIslands'

## dir.bombus <-
##     '/Volumes/bombus/Dropbox (University of Oregon)/skyIslands'

setwd(file.path(dir.bombus, "dataPrep"))

src.dir <- '../../skyIslands_saved/data/relational/relational/traditional/'
spec <-
    read.csv(file.path(src.dir, "specimens-complete.csv"),
             stringsAsFactors=FALSE)

source("src/misc.R")
source("src/prepNets.R")
source("src/specialization.R")


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

## subset to Net data only, pans have not been IDed from 2021+ as of
## Nov 2022
spec <- spec[spec$Method == "Net",]

## all possible sample round, site, combos of net data, used for veg
## data subsequently
collections <- data.frame(unique(cbind(spec$Site, spec$Year,
                                       spec$SampleRound)))
colnames(collections) <- c("Site", "Year", "SampleRound")

## drop non-bees, syrphids have been IDed up until 2021 (as of Nov
## 2022), butterflies, wasps need work

spec <- spec[spec$Family %in% c("Andrenidae", "Apidae",
                                 "Colletidae", "Halictidae",
                                "Megachilidae"),]

## ## drop non bees but keep the syrphids
## spec <- spec[spec$Family %in% c("Andrenidae", "Apidae",
##                                  "Colletidae", "Halictidae",
##                                  "Megachilidae", "Syrphidae"),]

## drop the the 2017 sample of PL because it was on fire for other
## sampling rounds and there was basically nothing blooming the first
## round, or leave it for phenology? As of Nov 2022 included because
## enough species have been IDed and useful as a early season data
## point
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

## check plant names
## library(Taxonstand)
## checked.plant.names <- TPL(id(spec$PlantGenusSpecies))

## write.csv(checked.plant.names,
## file="../../skyIslands_saved/data/checks/plant_names_check.csv")

spec.checked.plant.names <-
    read.csv(file="../../skyIslands_saved/data/checks/plant_names_check.csv")

spec <- fixPlantNames(spec, "PlantGenusSpecies", spec.checked.plant.names)

## ##  variabile identification fo Erigeron between years, combine
## to Erigerson spp.?
## spec$PlantGenus <- sapply(strsplit(spec$PlantGenusSpecies, " "),
##                           function(x) x[1])
## spec$PlantGenusSpecies[spec$PlantGenus == "Erigeron"] <-
##     "Erigeron spp. NA"

## *****************************************************************
## specimen-level parasite calculations
## *****************************************************************

dir.create("../data", showWarnings = FALSE)

parasites <- c( "AscosphaeraSpp",
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

save(spec, file="../data/spec_all_methods.Rdata")
write.csv(spec, file="../data/spec_all_methods.csv", row.names=FALSE)


## for networks, drop specimens withoutplant IDs, which will also drop
## pans and vanes
spec <- spec[spec$PlantGenusSpecies != "",]
spec <- spec[spec$GenusSpecies != "",]
save(spec, file="../data/spec_net.Rdata")
write.csv(spec, file="../data/spec_net.csv", row.names=FALSE)

## ***********************************************************************
## site/species level insect data
## ***********************************************************************

site.sp <- spec %>%
    group_by(Site, Year, SampleRound, GenusSpecies) %>%
    summarise(Abundance = length(GenusSpecies),
              SpParasitismRate=mean(ParasitePresence, na.rm=TRUE))

site.sum <- spec %>%
    group_by(Site, Year, SampleRound) %>%
    summarise(PollAbundance = length(GenusSpecies),
              PollRichness= length(unique(GenusSpecies)),
              VisitedFloralRichness= length(unique(PlantGenusSpecies)),
              BombusRichness= length(unique(GenusSpecies[Genus == "Bombus"])),
              PollDiversity=vegan:::diversity(table(GenusSpecies),
                                              index="shannon"),
              VisitedFloralDiversity=vegan:::diversity(table(PlantGenusSpecies),
                                              index="shannon"),
              BombusDiversity=vegan:::diversity(table(GenusSpecies[Genus == "Bombus"]),
                                                index="shannon"),
              SiteParasitismRate=mean(ParasitePresence, na.rm=TRUE),
              MeanParasiteRichness=mean(ParasiteRichness, na.rm=TRUE),
              SRDoyPoly1=mean(DoyPoly1),
              SRDoyPoly2=mean(DoyPoly2),
              SRDoy=mean(Doy),
              HBAbundance = sum(GenusSpecies == "Apis mellifera"),
              BombusAbundance = sum(Genus == "Bombus"),
              NonBombusHBAbundance =
                  sum(Genus != "Bombus" & Genus != "Apis"),
              HBSiteParasitismRate=mean(
                  ParasitePresence[GenusSpecies == "Apis mellifera"],
                  na.rm=TRUE),
              BombusSiteParasitismRate=mean(
                  ParasitePresence[Genus == "Bombus"], na.rm=TRUE))

## hb <- spec[spec$GenusSpecies == "Apis mellifera",]

## hb.site.sum <- hb %>%
##     group_by(Site, Year, SampleRound) %>%
##     summarise(HBAbundance =n(),
##               HBSiteParasitismRate=mean(ParasitePresence, na.rm=TRUE))

## site.sum  <- merge(site.sum, hb.site.sum, all.x=TRUE)

site.sp.yr <- spec %>%
    group_by(Site, Year, GenusSpecies, Genus) %>%
    summarise(Abundance = length(GenusSpecies))

bombus <- site.sp.yr[site.sp.yr$Genus == "Bombus",]
bombus$Genus  <- NULL

## add site characteristics
sites <- unique(data.frame(Site=spec$Site,
                           Lat= spec$Lat,
                           Area=spec$Area,
                           Elev=spec$Elev))
site.sum <- merge(site.sum, sites)
site.sum$Year <- as.factor(site.sum$Year)

## write species-level sumary data
write.csv(bombus, file='../data/bombus_year_site.csv', row.names=FALSE)
write.csv(site.sp.yr, file='../data/sp_year_site.csv', row.names=FALSE)
write.csv(site.sp, file='../data/spstats.csv', row.names=FALSE)

## write the site, year, sampling round summary data after merging
## with plant data

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

write.csv(traits, file='../data/networks_traits.csv', row.names=FALSE)

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
##  create site, SR, year level networks for plant-pollinators,
##  pollinator-parasites
## *******************************************************************

## plant-pollinator networks
makeNets(spec, net.type="YrSR")
makeNets(spec, net.type="Yr", mean.by.year=TRUE)

spec.sub <- agg.spec.sub %>%
    select(UniqueID, GenusSpecies, Site, Year, SampleRound,
        "AscosphaeraSpp",
           "ApicystisSpp", "CrithidiaExpoeki", "CrithidiaBombi",
            "NosemaBombi", "NosemaCeranae")

prep.para <- spec.sub %>%
    pivot_longer(cols=c("AscosphaeraSpp",
                        "ApicystisSpp", "CrithidiaExpoeki",
                        "CrithidiaBombi",
                         "NosemaBombi", "NosemaCeranae"),
                 names_to = "Parasite", values_to = "count")
prep.para <- as.data.frame(prep.para)

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
##  Data checks
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
##  Veg and bloom cleaning
## *******************************************************************

bloom <-
    read.csv(file.path(src.dir, "bloom-complete.csv"),
             stringsAsFactors=FALSE)
veg <-
    read.csv(file.path(src.dir, "veg-complete.csv"),
             stringsAsFactors=FALSE)

veg$PlantGenusSpecies <-  fix.white.space(paste(veg$PlantGenus,
                                          veg$PlantSpecies,
                                          veg$PlantVar))

## people in the field put NA as a placehold where there were no
## flowers to express data was collected, just nothing was
## there. Not needed for site level calculations.
veg <- veg[veg$PlantGenusSpecies != "",]

bloom$PlantGenusSpecies <-  fix.white.space(paste(bloom$PlantGenus,
                                          bloom$PlantSpecies,
                                          bloom$PlantVar))
bloom <- bloom[bloom$PlantGenusSpecies != "",]

## check plant names, run Nov 2022
## library(Taxonstand)
## bloom.checked.plant.names <- TPL(id(bloom$PlantGenusSpecies))
## write.csv(bloom.checked.plant.names,
## file="../../skyIslands_saved/data/checks/bloom_plant_names_check.csv")

## veg.checked.plant.names <- TPL(id(veg$PlantGenusSpecies))
## write.csv(veg.checked.plant.names,
## file="../../skyIslands_saved/data/checks/veg_plant_names_check.csv")

## update plant names
veg.checked.plant.names <-
    read.csv(file="../../skyIslands_saved/data/checks/veg_plant_names_check.csv")

bloom.checked.plant.names <-
    read.csv(file="../../skyIslands_saved/data/checks/bloom_plant_names_check.csv")

veg <- fixPlantNames(veg, "PlantGenusSpecies",
                     veg.checked.plant.names)
bloom <- fixPlantNames(bloom, "PlantGenusSpecies",
                     bloom.checked.plant.names)

id(veg$PlantGenusSpecies)
id(bloom$PlantGenusSpecies)

## fix dates
veg$Date <- as.Date(veg$Date,  "%Y-%m-%d")
veg$Year <- format(veg$Date,  "%Y")

bloom$Date <- as.Date(bloom$Date,  "%Y-%m-%d")
bloom$Year <- format(bloom$Date,  "%Y")

## 2017-2022 quad options
quads.2017.2022 <- unique(veg$Quadrat)
quads.2012 <- unique(veg$Quadrat[veg$Year == 2012])
## quads not done in 2012
not.in.2012 <- quads.2017.2022[!quads.2017.2022 %in% quads.2012]

## prep a matrix of all the quads for site richness etc.
combos <- collections[rep(seq_len(nrow(collections)),
                          each = length(quads.2017.2022)), ]
combos$Quadrat <- quads.2017.2022
## drop the quads not done in 2012
combos <- combos[!(combos$Year == 2012 & combos$Quad %in%
                   not.in.2012),]

## blooming flowers
veg.blooming <- veg[veg$NumBlooms != "" &
                    veg$NumBlooms != "0",]

## didn't count the exact number of flowers in 2017-2018, use
## midpoints of bins
veg.blooming$NumBloomsNum <- veg.blooming$NumBlooms
veg.blooming$NumBloomsNum[veg.blooming$NumBlooms == "<10"] <- 5
veg.blooming$NumBloomsNum[veg.blooming$NumBlooms == "10-100"] <- 50
veg.blooming$NumBloomsNum[veg.blooming$NumBlooms == "100-1000"] <- 500

veg.blooming$NumBloomsNum <- as.numeric(veg.blooming$NumBloomsNum)

## *******************************************************************
##  veg site summaries as mean of quads
## *******************************************************************

## by quadrat for a site average
veg.bloom.quad.sp <- veg.blooming %>%
    group_by(Quadrat, Site, Year, PlantGenusSpecies, SampleRound) %>%
    summarise(FloralAbundance = sum(NumBloomsNum),
              FloweringPlantAbundance= sum(PlantCount, na.rm=TRUE))

veg.bloom.quad.sp$KeyQuad <- paste(veg.bloom.quad.sp$Site,
                              veg.bloom.quad.sp$Quadrat,
                              veg.bloom.quad.sp$Year,
                              veg.bloom.quad.sp$SampleRound)

veg.div.quad <- tapply(veg.bloom.quad.sp$FloralAbundance,
                  veg.bloom.quad.sp$KeyQuad,
                  vegan::diversity)

veg.plant.div.quad <- tapply(veg.bloom.quad.sp$FloweringPlantAbundance,
                  veg.bloom.quad.sp$KeyQuad,
                  vegan::diversity)


veg.bloom.sum.quad <- veg.bloom.quad.sp %>%
    group_by(Quadrat, Site, Year, SampleRound) %>%
    summarise(FloralAbundance = sum(FloralAbundance),
              FloralRichness= length(unique(PlantGenusSpecies)),
              FloweringPlantAbundance = sum(FloweringPlantAbundance)
              )

veg.name <-   names(veg.div.quad)
names(veg.div.quad) <- NULL

## add div data into site summary data
veg.bloom.sum.quad$FloralDiversity <-  veg.div.quad[match(
    paste(veg.bloom.sum.quad$Site,
          veg.bloom.sum.quad$Quadrat,
          veg.bloom.sum.quad$Year,
          veg.bloom.sum.quad$SampleRound),
    veg.name)]

veg.bloom.sum.quad$FloweringPlantDiversity <-
    veg.plant.div.quad[match(
        paste(veg.bloom.sum.quad$Site,
              veg.bloom.sum.quad$Quadrat,
              veg.bloom.sum.quad$Year,
              veg.bloom.sum.quad$SampleRound),
        veg.name)]

## merge together quad combos and quad level data
combos <- merge(combos, veg.bloom.sum.quad, all.x=TRUE)

cols.to.fill <-  c("FloralAbundance", "FloralRichness",
                   "FloralDiversity", "FloweringPlantAbundance",
                   "FloweringPlantDiversity")
## set NAs to zero
combos[, cols.to.fill][is.na(combos[, cols.to.fill])] <- 0

veg.bloom.mean <- combos %>%
    group_by(Site, Year, SampleRound) %>%
    summarise(MeanFloralAbundance = mean(FloralAbundance),
              MeanFloralRichness= mean(FloralRichness),
              MeanFloralDiversity =mean(FloralDiversity),
              MeanFloweringPlantDiversity =mean(FloweringPlantDiversity),
              MeanFloweringPlantAbundance =mean(FloweringPlantAbundance)
              )

veg.bloom.mean[order(veg.bloom.mean$MeanFloralRichness), ]

write.csv(veg.bloom.mean, file="../data/veg.csv", row.names=FALSE)

## merge with specimen summary data
site.sum <- merge(site.sum, veg.bloom.mean, all.x=TRUE)
write.csv(site.sum, file='../data/sitestats.csv', row.names=FALSE)

## floral richness across the entire meadow across sampling rounds
veg.year.sum <- veg.blooming %>%
    group_by(Site, Year) %>%
    summarise(TotalFloralRichness= length(unique(PlantGenusSpecies)))

write.csv(veg.year.sum, file="../data/veg_species_richness.csv", row.names=FALSE)

## not using site-level bloom data currently, has not been checked
## over as of Nov 2022

## *******************************************************************
## checking veg data between years
## *******************************************************************

## to.sample <- c("CH", "HM", "JC", "MM", "PL", "RP", "SC", "SM")

## veg.prep <- veg[veg$Site %in% to.sample,]
## veg.split <- split(veg.prep, veg.prep$Site)

## spec.veg <- spec[spec$Site %in% to.sample,]
## spec.split <- split(spec.veg, spec.veg$Site)

## checkveg <- function(x){
##     print(unique(x$Site))
##     years <- tapply(x$PlantGenusSpecies, x$Year, unique)
##     if(length(years) == 2){
##         badmatch12 <- years[[1]][!years[[1]] %in% years[[2]]]
##         badmatch21 <- years[[2]][!years[[2]] %in% years[[1]]]
##         return(list(badmatch12=badmatch12,
##                     badmatch21=badmatch21))
##     } else if(length(years) == 3){
##         badmatch12 <- years[[1]][!years[[1]] %in% years[[2]]]
##         badmatch21 <- years[[2]][!years[[2]] %in% years[[1]]]

##         badmatch13 <- years[[1]][!years[[1]] %in% years[[3]]]
##         badmatch31 <- years[[3]][!years[[3]] %in% years[[1]]]

##         badmatch23 <- years[[2]][!years[[2]] %in% years[[3]]]
##         badmatch32 <- years[[3]][!years[[3]] %in% years[[2]]]

##         return(list(badmatch12=badmatch12,
##                     badmatch21=badmatch21,
##                     badmatch13= badmatch13,
##                     badmatch31=badmatch31,
##                     badmatch23=badmatch23,
##                     badmatch32=badmatch32))
##     } else{
##         return(NULL)
##     }
## }

## badmatches.veg <- lapply(veg.split, checkveg)

## badmatches.spec <- lapply(spec.split, checkveg)


## all.plant.sp <- rbind(data.frame(Site=veg.prep$Site,
##                                  PlantGenusSpecies=
##                                      veg.prep$PlantGenusSpecies),
##                       data.frame(Site=spec.veg$Site,
##                                  PlantGenusSpecies=
##                                      spec.veg$PlantGenusSpecies))

## all.plant.sp <- unique(all.plant.sp)
## all.plant.sp <- sort(all.plant.sp)

## write.csv(all.plant.sp, row.names=FALSE,
##           file="../data/vegbySite.csv")
