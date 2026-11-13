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
setwd("~/")
# setwd('C:/')
source("lab_paths.R")
local.path

dir.bombus <- file.path(local.path, "skyIslands")

## *****************************************************************
## create relational database, add species IDs
## *****************************************************************

setwd(dir.bombus)
source('dataPrep/relational/1prep.R')

setwd(dir.bombus)
source('dataPrep/relational/1_2prep_climate.R')

setwd(dir.bombus)
source('dataPrep/relational/2make.R')

setwd(dir.bombus)
source('dataPrep/relational/3join.R')

## *****************************************************************
## prep specimen data
## *****************************************************************

rm(list=ls())
setwd("~/")
# setwd("C:/")
source("lab_paths.R")
setwd(file.path(local.path, "skyIslands_saved"))

spec <-
    read.csv('data/relational/relational/traditional/specimens-complete.csv',
             stringsAsFactors=FALSE)

setwd("../skyIslands/dataPrep")
source("src/misc.R")
source("src/join_checks.R")
source("src/prepNets.R")
source("src/guild_networks.R")
source("src/specialization.R")
source("src/SpCalcMetrics.R")
source("src/vaznull2.R")

spec$GenusSpecies <- fix.white.space(paste(
  ifelse(is.na(spec$Genus), '', spec$Genus),
  ifelse(is.na(spec$Species), '', spec$Species),
  ifelse(is.na(spec$SubSpecies), '', spec$SubSpecies)))

if (!'PlantGenusSpecies' %in% names(spec)) {
  spec$PlantGenusSpecies <- fix.white.space(paste(
    ifelse(is.na(spec$PlantGenus), '', spec$PlantGenus),
    ifelse(is.na(spec$PlantSpecies), '', spec$PlantSpecies),
    ifelse(is.na(spec$PlantVar), '', spec$PlantVar)))
}

missing.plant <- is.na(spec$PlantGenusSpecies) |
                 spec$PlantGenusSpecies == ''
if ('SourcePlantName' %in% names(spec)) {
  spec$PlantGenusSpecies[missing.plant] <-
    spec$SourcePlantName[missing.plant]
}

spec$Date <- as.Date(spec$Date, format='%Y-%m-%d')
spec$Doy <- as.numeric(strftime(spec$Date, format='%j'))
spec$Year <- as.numeric(format(spec$Date,'%Y'))

## all possible sample round, site, combos of net data, used for veg
## data subsequently
collections <- data.frame(unique(cbind(spec$Site, spec$Year,
                                       spec$SampleRound)))
colnames(collections) <- c("Site", "Year", "SampleRound")

## All plant names already come from the reconciled relational plant
## table; no second taxonomy pass is applied during specimen prep.

spec$Int <- fix.white.space(paste(spec$GenusSpecies,
                                  spec$PlantGenusSpecies))
check_unique_key(spec, 'UniqueID', 'Specimens loaded for analysis')

## *****************************************************************
## specimen-level parasite calculations
## *****************************************************************
dir.create("../data/", showWarnings = FALSE)
dir.create("../data/networks", showWarnings = FALSE)
dir.create("../data/splevel_network_metrics", showWarnings = FALSE)

crithidias <- c("CrithidiaExpoeki",
                "CrithidiaBombi", "CrithidiaSpp",
                "CrithidiaMellificae")
nosemas <- c("NosemaBombi", "NosemaCeranae")
parasites <- c( "AscosphaeraSpp",
               "ApicystisSpp", crithidias, nosemas)

## Apidae == 1 means successful DNA extraction.
## Only these specimens can be coded 0 for parasite tests with no positive
## call. Their 0 values are established upstream in 1prep.R. Explicit
## NA flags for contaminated/uninterpretable assays MUST remain NA.
## Apidae != 1 or missing is unscreened: no parasite-negative zeros.
for (this.para in parasites) {
  spec[[this.para]] <- suppressWarnings(as.numeric(spec[[this.para]]))
}
spec$Apidae <- suppressWarnings(as.numeric(spec$Apidae))
not.screened <- is.na(spec$Apidae) | spec$Apidae != 1
spec[not.screened, parasites] <- NA
spec$DNAExtractionSuccessful <- ifelse(spec$Apidae == 1, 1L, NA_integer_)
##  If a successful Apidae assay still has NA parasite results,
## retain NA: these include explicit contaminated-assay exclusions.
## Do not turn an indeterminate assay into a documented negative.

## Denominator for pathogen-specific prevalence is the number of
## Partial assay panels have an indeterminate *total* richness.
known.richness <- function(x) {
  if (anyNA(x)) return(NA_real_)
  sum(x)
}
## A positive is certain even if another assay is missing; a negative is
## certain only if every assay in this set is a documented negative.
known.presence <- function(x) {
  if (any(x > 0, na.rm=TRUE)) return(1L)
  if (anyNA(x)) return(NA_integer_)
  0L
}
spec$ParasiteRichness <- apply(spec[, parasites, drop=FALSE], 1,
                               known.richness)
spec$CrithidiaRichness <- apply(spec[, crithidias, drop=FALSE], 1,
                                known.richness)
spec$NosemaRichness <- apply(spec[, nosemas, drop=FALSE], 1,
                             known.richness)
spec$PossibleParasite <- rowSums(!is.na(spec[, parasites, drop=FALSE]))
spec$ParasitePresence <- apply(spec[, parasites, drop=FALSE], 1,
                               known.presence)
spec$CrithidiaPresence <- apply(spec[, crithidias, drop=FALSE], 1,
                                known.presence)
spec$NosemaPresence <- apply(spec[, nosemas, drop=FALSE], 1,
                             known.presence)
check.spec <- spec[!is.na(spec$Apidae),]
check.spec <- check.spec[check.spec$GenusSpecies == "",]

write.csv(check.spec,
 file="../../skyIslands_saved/data/checks/screened_no_ID.csv",
 row.names=FALSE)


write.csv(spec, file="../data/spec_all_methods.csv", row.names=FALSE)

spec.net <- spec[spec$Method == "Net",]
spec.pan <- spec[spec$Method == "Pan",]
spec <- spec[spec$Method != "Vane",]

write.csv(spec, file="../data/spec_net_pan.csv", row.names=FALSE)

## ***********************************************************************
## site/species level insect data
## ***********************************************************************

bee.families <- c("Andrenidae", "Apidae", "Colletidae", "Halictidae",
                  "Megachilidae", "Melittidae")

## Butterfly families plus Sphingidae (hawk moths).
## Other moth families (e.g., Erebidae, Noctuidae) are not included.

butterfly.hawk.families <- c('Hesperiidae', 'Lycaenidae', 'Nymphalidae',
                            'Papilionidae', 'Pieridae', 'Riodinidae',
                            'Sphingidae')
focal.pollinator.families <- c(bee.families, 'Syrphidae',
                               butterfly.hawk.families)
## Keep the existing Bees+Syrphids+Lep columns unchanged. Add a separate
## combined guild that includes bee flies (Bombyliidae).
all.pollinator.families <- c(focal.pollinator.families, 'Bombyliidae')


calcSummaryStats <- function(spec.method, method){
    site.sp <- spec.method %>%
        group_by(Site, Year, SampleRound, GenusSpecies) %>%
        summarise(Abundance = length(GenusSpecies),
                  SpParasitismRate = mean(ParasitePresence, na.rm=TRUE),
                  SpCrithidiaParasitismRate = mean(CrithidiaPresence, na.rm=TRUE),
                  SpApicystisParasitismRate = mean(ApicystisSpp, na.rm=TRUE),
                  SpCrithidiaBombusParasitismRate = mean(CrithidiaPresence[Genus == "Bombus"],
                                                      na.rm=TRUE),
                  SpCrithidiaHBParasitismRate = mean (CrithidiaPresence
                                                   [GenusSpecies == "Apis mellifera"],
                                                   na.rm=TRUE),
                  SpParasitism = sum(ParasitePresence, na.rm=TRUE),
                  SpCrithidiaPresence= sum(CrithidiaPresence, na.rm=TRUE),
                  SpApicystisSpp = sum(ApicystisSpp, na.rm=TRUE),
                  SpNosemaBombi= sum(NosemaBombi, na.rm=TRUE),
                  SpNosemaCeranae = sum(NosemaCeranae, na.rm=TRUE),
                  SpAscosphaeraSpp= sum(AscosphaeraSpp, na.rm=TRUE),
                  SpCrithidiaExpoeki = sum(CrithidiaExpoeki, na.rm=TRUE),
                  SpCrithidiaMellificae = sum(CrithidiaMellificae, na.rm=TRUE),
                  SpCrithidiaBombi = sum(CrithidiaBombi, na.rm=TRUE),
                  SpCrithidiaSpp = sum(CrithidiaSpp, na.rm=TRUE),
                  SpScreened = sum(Apidae == 1, na.rm=TRUE)
                  )

    ## Preserve all existing summary-statistic
    ## definitions/columns and add new guilds below. 
    ## Syrphidae, and butterflies + Sphingidae but not other fly families.
    site.sum <- spec.method %>%
        group_by(Site, Year, SampleRound) %>%
        summarise(PollAbundance = length(GenusSpecies),
                  BeeAbundance = length(GenusSpecies[Family %in%
                                                     bee.families]),
                  SyrphidAbundance =
                      length(GenusSpecies[Family == 'Syrphidae']),
                  BeeFlyAbundance = sum(Family %in% 'Bombyliidae'),
                  HBAbundance = sum(GenusSpecies == 'Apis mellifera'),
                  BombusAbundance = sum(Genus == 'Bombus'),
                  NonBombusHBAbundance =
                      sum(Genus != 'Bombus' &
                          Genus != 'Apis' &
                          Family %in% bee.families),
                  ## New Lepidoptera and combined counts.
                  ButterflyHawkAbundance =
                      sum(Family %in% butterfly.hawk.families),
                  BeeSyrphLepAbundance =
                      sum(Family %in% focal.pollinator.families),
                  BeeSyrphBeeFlyLepAbundance =
                      sum(Family %in% all.pollinator.families),

                  PollRichness= length(unique(GenusSpecies)),
                  BeeRichness= length(unique(GenusSpecies[Family %in%
                                                          bee.families])),
                  SyrphidRichness= length(unique(
                    GenusSpecies[Family  == 'Syrphidae'])),
                  BeeFlyRichness = guild.richness(
                    GenusSpecies[Family %in% 'Bombyliidae']),
                  BombusRichness= length(unique(
                    GenusSpecies[Genus == 'Bombus'])),
                  ## Lepidoptera and
                  ## combined taxonomic-label richness. 
                  ButterflyHawkRichness = guild.richness(
                      GenusSpecies[Family %in% butterfly.hawk.families]),
                  BeeSyrphLepRichness = guild.richness(
                      GenusSpecies[Family %in% focal.pollinator.families]),
                  BeeSyrphBeeFlyLepRichness = guild.richness(
                      GenusSpecies[Family %in% all.pollinator.families]),

                  PollDiversity=vegan::diversity(table(GenusSpecies),
                                                  index='shannon'),
                  BeeDiversity=vegan::diversity(table(
                    GenusSpecies[Family %in% bee.families]),
                                                 index='shannon'),
                  SyrphidDiversity=vegan::diversity(table(
                    GenusSpecies[Family  == 'Syrphidae']),
                                                     index='shannon'),
                  BeeFlyDiversity = guild.shannon(
                    GenusSpecies[Family %in% 'Bombyliidae']),
                  BombusDiversity=vegan::diversity(table(
                    GenusSpecies[Genus == 'Bombus']),
                                                    index='shannon'),
                  ## Shannon for three added groups.
                  ButterflyHawkDiversity = guild.shannon(
                      GenusSpecies[Family %in% butterfly.hawk.families]),
                  BeeSyrphLepDiversity = guild.shannon(
                      GenusSpecies[Family %in% focal.pollinator.families]),
                  BeeSyrphBeeFlyLepDiversity = guild.shannon(
                      GenusSpecies[Family %in% all.pollinator.families]),

                  VisitedFloralRichness= length(unique(PlantGenusSpecies)),
                  VisitedFloralDiversity=vegan::diversity(table(
                    PlantGenusSpecies),  index='shannon'),

                  SiteParasitismRate=mean(ParasitePresence, na.rm=TRUE),
                  MeanParasiteRichness=mean(ParasiteRichness, na.rm=TRUE),
                  HBSiteParasitismRate=mean(
                      ParasitePresence[GenusSpecies == 'Apis mellifera'],
                      na.rm=TRUE),
                  BombusSiteParasitismRate=mean(
                      ParasitePresence[Genus == 'Bombus'],
                    na.rm=TRUE),
                  SRDoy=mean(Doy),
                  CrithidiaParasitismRate=mean(CrithidiaPresence, na.rm=TRUE),
                  ApicystisParasitismRate=mean(ApicystisSpp, na.rm=TRUE),
                  CrithidiaBombusParasitismRate= mean(CrithidiaPresence[Genus == 'Bombus'],
                                                      na.rm=TRUE),
                  CrithidiaHBParasitismRate= mean (CrithidiaPresence
                                                   [GenusSpecies == 'Apis mellifera'],
                                                   na.rm=TRUE))

    site.sp.yr <- spec.method %>%
        group_by(Site, Year, GenusSpecies, Genus) %>%
        summarise(Abundance = length(GenusSpecies))


    site.sp.yr.round <- spec.method %>%
      group_by(Site, Year, SampleRound, GenusSpecies, Genus) %>%

      summarise(AbundanceSYR = length(GenusSpecies))

    bombus <- site.sp.yr[site.sp.yr$Genus == "Bombus",]
    bombus$Genus  <- NULL

    ## write species-level sumary data
    write.csv(bombus, file=sprintf('../data/bombus_year_site_%s.csv',
                                   method),
              row.names=FALSE)
    write.csv(site.sp.yr, file=sprintf('../data/sp_year_site_%s.csv',
                                       method),
              row.names=FALSE)
    write.csv(site.sp, file=sprintf('../data/spstats_%s.csv', method),
              row.names=FALSE)

    write.csv(site.sp.yr.round,
              file=sprintf('../data/sp_year_site_round_%s.csv', method),
              row.names=FALSE)
    if (method == 'all') write.csv(site.sp.yr.round,
      file='../data/sp_year_site_round.csv', row.names=FALSE)
    return(site.sum)
}

net.site.sum <- calcSummaryStats(spec.net, "net")
pan.site.sum <- calcSummaryStats(spec.pan, "pan")
all.site.sum <- calcSummaryStats(spec, "all")

## ***********************************************************************
## NETTED genus-by-year and ID-completeness diagnostics.
## Export diagnostic CSVs and plots to skyIslands_saved/data/checks,
##  Source is netted specimens before any
## summary joins, so each row is counted once.
## ***********************************************************************
dir.guild.audit <- '../../skyIslands_saved/data/checks/guild_id_review'
dir.create(dir.guild.audit, recursive=TRUE, showWarnings=FALSE)
source("src/guide_review.R")

table(all.site.sum$Site)
table(net.site.sum$Site)
table(pan.site.sum$Site)

##  Include newly calculated guild columns in the
## existing summary-variable reference list.
sum.cols <- c("PollAbundance", "BeeAbundance", "SyrphidAbundance",
              "BeeFlyAbundance", "BeeFlyRichness", "BeeFlyDiversity",
              "BeeSyrphBeeFlyLepAbundance", "BeeSyrphBeeFlyLepRichness",
              "BeeSyrphBeeFlyLepDiversity",
              "HBAbundance", "BombusAbundance", "ButterflyHawkAbundance",
              "BeeSyrphLepAbundance", "NonBombusHBAbundance",
              "PollRichness", "BeeRichness", "SyrphidRichness",
              "BombusRichness", "HBRichness", "ButterflyHawkRichness",
              "BeeSyrphLepRichness", "PollDiversity", "BeeDiversity",
              "SyrphidDiversity", "BombusDiversity", "HBDiversity",
              "ButterflyHawkDiversity", "BeeSyrphLepDiversity",
              "VisitedFloralRichness", "VisitedFloralDiversity",
              "SiteParasitismRate", "MeanParasiteRichness",
              "HBSiteParasitismRate", "BombusSiteParasitismRate")


keys.sr <- c('Site', 'Year', 'SampleRound')
dup.cols <- c('SRDoyPoly1', 'SRDoyPoly2', 'SRDoy')
net.site.sum <- net.site.sum[, !names(net.site.sum) %in% dup.cols]
pan.site.sum <- pan.site.sum[, !names(pan.site.sum) %in% dup.cols]
net.vars <- !names(net.site.sum) %in% keys.sr
pan.vars <- !names(pan.site.sum) %in% keys.sr
names(net.site.sum)[net.vars] <- paste0('Net_', names(net.site.sum)[net.vars])
names(pan.site.sum)[pan.vars] <- paste0('Pan_', names(pan.site.sum)[pan.vars])
check_unique_key(all.site.sum, keys.sr, 'All-method sample-round stats')
check_unique_key(net.site.sum, keys.sr, 'Net sample-round stats')
check_unique_key(pan.site.sum, keys.sr, 'Pan sample-round stats')
site.sum <- dplyr::left_join(all.site.sum, net.site.sum, by=keys.sr,
                             relationship='one-to-one')
site.sum <- dplyr::left_join(site.sum, pan.site.sum, by=keys.sr,
                             relationship='one-to-one')
check_rows_preserved(all.site.sum, site.sum,
                     'Combined site x year x round stats')


## add back site characteristics
## Event-derived geography below fills values for empty sampling rounds.
mean.unique <- function(x) {
  x <- unique(suppressWarnings(as.numeric(x)))
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}
sites <- spec %>%
  group_by(Site, Year) %>%
  summarise(Lat=mean.unique(Lat), Area=mean.unique(Area),
            Elev=mean.unique(Elev), .groups='drop')
check_unique_key(sites, c('Site', 'Year'), 'Site-year metadata')
site.sum <- dplyr::left_join(site.sum, sites, by=c('Site', 'Year'),
                              relationship='many-to-one')
site.sum$Year <- as.factor(site.sum$Year)

## write the site, year, sampling round summary data after merging
## with plant data


## ***********************************************************************
## For multi-study parasite summaries, a bit redundant but nice to
## have in a standardized format
## ***********************************************************************

par.site.sp <- spec.net %>%
    group_by(Site, Year, SampleRound, GenusSpecies, Genus) %>%
    summarise(SpAbundance = length(GenusSpecies),
              SpParasitism = sum(ParasitePresence, na.rm=TRUE),
              SpCrithidiaPresence= sum(CrithidiaPresence, na.rm=TRUE),
              SpApicystisSpp = sum(ApicystisSpp, na.rm=TRUE),
              SpNosemaBombi= sum(NosemaBombi, na.rm=TRUE),
              SpNosemaCeranae = sum(NosemaCeranae, na.rm=TRUE),
              SpAscosphaeraSpp= sum(AscosphaeraSpp, na.rm=TRUE),
              SpCrithidiaExpoeki = sum(CrithidiaExpoeki, na.rm=TRUE),
              SpCrithidiaBombi = sum(CrithidiaBombi, na.rm=TRUE),
              SpCrithidiaSpp = sum(CrithidiaSpp, na.rm=TRUE),
              SpScreened = sum(Apidae == 1, na.rm=TRUE)
              )

par.site.gen <- spec.net %>%
    group_by(Site, Year, SampleRound, Genus) %>%
    summarise(GenusAbundance = length(GenusSpecies),
              GenusParasitism = sum(ParasitePresence, na.rm=TRUE),
              GenusCrithidiaPresence= sum(CrithidiaPresence, na.rm=TRUE),
              GenusApicystisSpp = sum(ApicystisSpp, na.rm=TRUE),
              GenusNosemaBombi= sum(NosemaBombi, na.rm=TRUE),
              GenusNosemaCeranae = sum(NosemaCeranae, na.rm=TRUE),
              GenusAscosphaeraSpp= sum(AscosphaeraSpp, na.rm=TRUE),
              GenusCrithidiaExpoeki = sum(CrithidiaExpoeki, na.rm=TRUE),
              GenusCrithidiaBombi = sum(CrithidiaBombi, na.rm=TRUE),
              GenusCrithidiaSpp = sum(CrithidiaSpp, na.rm=TRUE),
              GenusScreened = sum(Apidae == 1, na.rm=TRUE)
              )

par.site.sum <- spec.net %>%
    group_by(Site, Year, SampleRound) %>%
    summarise(SiteBeeAbundance = length(GenusSpecies[Family %in%
                                                     bee.families]),

              SiteBeeRichness= length(unique(GenusSpecies[Family %in%
                                                          bee.families])),

              SiteBeeDiversity=vegan::diversity(table(
                                           GenusSpecies[Family %in% bee.families]),
                                           index="shannon"),

              SiteParasitism = sum(ParasitePresence, na.rm=TRUE),
              SiteCrithidiaPresence= sum(CrithidiaPresence, na.rm=TRUE),
              SiteApicystisSpp = sum(ApicystisSpp, na.rm=TRUE),
              SiteNosemaBombi= sum(NosemaBombi, na.rm=TRUE),
              SiteNosemaCeranae = sum(NosemaCeranae, na.rm=TRUE),
              SiteAscosphaeraSpp= sum(AscosphaeraSpp, na.rm=TRUE),
              SiteCrithidiaExpoeki = sum(CrithidiaExpoeki, na.rm=TRUE),
              SiteCrithidiaBombi = sum(CrithidiaBombi, na.rm=TRUE),
              SiteCrithidiaSpp = sum(CrithidiaSpp, na.rm=TRUE),

              SiteScreened = sum(Apidae == 1, na.rm=TRUE)
              )



all.sums <- left_join(par.site.sp, par.site.gen,
  by=c('Site', 'Year', 'SampleRound', 'Genus'),
  relationship='many-to-one')
all.sums <- left_join(all.sums, par.site.sum,
  by=c('Site', 'Year', 'SampleRound'),
  relationship='many-to-one')
all.sums$Project <- "SI"

save(all.sums, file="../data/SI_parasite_sums.Rdata")


## *******************************************************************
## create a giant plant-pollinator network to calculate specialization
## etc. across all SI
## ***************************************** **************************
##  exclude empty plant or pollinator IDs from observed interactions.
## Apply the same basic ID consistency checks before calculating
## dataset-wide plant-pollinator traits and event-level guild networks.
## Raw specimen records themselves are never deleted.
spec.net.nets <- audit_network_ids(spec.net)$usable

agg.spec <- aggregate(list(abund=spec.net.nets$GenusSpecies),
                      list(GenusSpecies=spec.net.nets$GenusSpecies,
                           PlantGenusSpecies=spec.net.nets$PlantGenusSpecies),
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
agg.spec.sub <- spec[spec$Apidae == 1 & !is.na(spec$Apidae),]

agg.spec.para <- aggregate(agg.spec.sub[, parasites],
                      list(GenusSpecies=agg.spec.sub$GenusSpecies),
                                              sum, na.rm=TRUE)

para.gensp.counts <- table(agg.spec.sub$GenusSpecies)
para.gensp.counts

para.gen.counts <- table(agg.spec.sub$Genus)
para.gen.counts

## proportion of individuals screened
agg.spec.para[, parasites] <- agg.spec.para[, parasites]/
    para.gensp.counts[agg.spec.para$GenusSpecies]

nets.para <- agg.spec.para
rownames(nets.para) <- nets.para$GenusSpecies
nets.para$GenusSpecies <- NULL
nets.para <- as.matrix(nets.para)


all.traits.para <- specieslevel(nets.para)
## calculate rarified plant.pol degree

para.traits <- data.frame(GenusSpecies= unlist(sapply(all.traits.para,
                                                 rownames)),
                     do.call(rbind, all.traits.para))

rownames(para.traits) <- NULL

write.csv(para.traits, file='../data/parasitetraits.csv', row.names=FALSE)

## *******************************************************************
##  create site, SR, year level networks for plant-pollinators,
##  pollinator-parasites
## *******************************************************************

## Construct an explicitly filtered named list of plant x insect matrices
## for every pollinator guild, at Site x Year x SampleRound only.
## The event export also records true netting rounds with zero observations.
events.for.nets <- read.csv(paste0(
  '../../skyIslands_saved/data/relational/relational/traditional/',
  'conditions-specimens-complete.csv'), stringsAsFactors=FALSE,
  colClasses=c(SampleRound='character'))
net.events <- events.for.nets[!is.na(events.for.nets$Method) &
  tolower(trimws(events.for.nets$Method)) == 'net',
  c('Site','Year','SampleRound'), drop=FALSE]
if (!nrow(net.events)) warning('No documented Net events: ',
  'networks will use observed specimen events only.')
## Reload the network helper when this block is run interactively after
## editing the file; otherwise a previous definition remains in memory.
source('src/guild_networks.R')
net.guild.result <- make_guild_networks(
  spec.net, net.events,
  guilds=pollinator_guilds(bee.families, butterfly.hawk.families),
  network_dir='../data/networks',
  metrics_dir='../data/splevel_network_metrics',
  audit_dir='../../skyIslands_saved/data/checks/network_id_review')
bees.yr.sr <- net.guild.result$bee_metrics

## ## Calculates species level network metrics by site and round
## ## Takes a while, only re-run if bee data is changed 
## N <- 99
## load('../data/networks/YearSR_PlantPollinator_Bees.Rdata')

## sp.mets <- lapply(nets, SpCalcNetworkMetrics,
##                N=N, index=c("closeness", "betweenness",
##                             "degree", "d", "normalised degree"))
## save(sp.mets,
##      file="../data/sp_network_mets_sr.RData")

#sp.mets <- sp.mets[!sapply(sp.mets, function(x) is.null(dim(x)))]
#sp.network.metrics <- SpPrepDat(sp.mets,  spec.net.nets,
#                   cols.to.keep=c("GenusSpecies", "Site","Year", "SampleRound"),
#                   net.type="YearSR")

## Parasites
spec.sub <- agg.spec.sub %>%
  dplyr::select(UniqueID, GenusSpecies, Site, Year, SampleRound, AscosphaeraSpp,
         ApicystisSpp, CrithidiaExpoeki, CrithidiaBombi,
         CrithidiaMellificae,
         NosemaBombi, NosemaCeranae)

prep.para <- spec.sub %>%
  pivot_longer(cols=c("AscosphaeraSpp",
                      "ApicystisSpp", "CrithidiaExpoeki",
                      "CrithidiaBombi", "CrithidiaMellificae",
                      "NosemaBombi", "NosemaCeranae"),
               names_to = "Parasite", values_to = "count")
prep.para <- as.data.frame(prep.para)

## an NA assay is unobserved, not an absent parasite edge;
## omit those particular assays before aggregating species prevalence.
prep.para <- prep.para[!is.na(prep.para$count), ]

par.bees.yr.sr <- makeNets(prep.para, net.type="YearSR",
                           species=c("Pollinator",
                                     "Parasite"),
                           lower.level="GenusSpecies",
                           higher.level="Parasite",
                           poll.groups="Bees")


## Annual mean parasite networks are intentionally not generated.
## Each sampling round's network is retained in par.bees.yr.sr above.

## merge site summary metrics
print("specimen dim before network metric year-sr merge")
dim(spec.net)

bees.yr.sr <- bees.yr.sr[bees.yr.sr$speciesType == 'higher.level',]
bees.yr.sr$SpSiteYear <- NULL
network.keys <- c('Site', 'Year', 'SampleRound', 'GenusSpecies')
bees.yr.sr$Year <- as.character(bees.yr.sr$Year)
spec.net$Year <- as.character(spec.net$Year)
bees.yr.sr$SampleRound <- as.character(bees.yr.sr$SampleRound)
spec.net$SampleRound <- as.character(spec.net$SampleRound)
check_unique_key(bees.yr.sr, network.keys, 'Bee network metrics')
spec.net.before <- spec.net
spec.net <- dplyr::left_join(spec.net, bees.yr.sr, by=network.keys,
                             relationship='many-to-one')
check_rows_preserved(spec.net.before, spec.net, 'Network metrics -> specimens')
dim(spec.net)
print("specimen dim after network metric year-sr merge")

spec.net$speciesType <- NULL

## *******************************************************************
##  Data checks
## *******************************************************************

print(paste("Pollinator species", length(unique(spec$GenusSpecies))))
print(paste("Bee species", length(unique(spec$GenusSpecies[
                                                spec$Family %in%
                                                bee.families]))))
print(paste("Syrphid species", length(unique(spec$GenusSpecies[
                                                spec$Family == "Syrphidae"]))))
print(paste("Plant species", length(unique(spec$PlantGenusSpecies))))
print(paste("Pollinator genera", length(unique(spec$Genus))))
print(paste("Interactions", length(unique(spec$Int))))
print(paste("Specimens", nrow(spec)))

## *******************************************************************
##  Veg and bloom cleaning
## *******************************************************************

bloom <-
    read.csv("../../skyIslands_saved/data/relational/relational/traditional/bloom-complete.csv",
             stringsAsFactors=FALSE)
veg <-
    read.csv("../../skyIslands_saved/data/relational/relational/traditional/veg-complete.csv",
             stringsAsFactors=FALSE)

veg.required <- c('Site', 'Date', 'SampleRound', 'Quadrat',
                  'PlantGenus', 'PlantSpecies', 'PlantVar',
                  'NumBlooms', 'PlantCount')
missing.veg <- setdiff(veg.required, names(veg))
if (length(missing.veg)) {
  stop('veg-complete.csv is missing: ', paste(missing.veg, collapse=', '),
       '. Available columns: ', paste(names(veg), collapse=', '),
       '. Check the current 2make.R / 3join.R vegetation export.')
}

veg$Date <- as.Date(veg$Date, format='%Y-%m-%d')
if (anyNA(veg$Date)) {
  stop('veg-complete.csv contains ', sum(is.na(veg$Date)),
       ' missing/unparseable vegetation dates; inspect its Date column.')
}
veg$Year <- as.numeric(format(veg$Date, '%Y'))
veg$SampleRound <- as.character(veg$SampleRound)
veg$Quadrat <- trimws(as.character(veg$Quadrat))

bloom$Date <- as.Date(bloom$Date, format='%Y-%m-%d')
bloom$Year <- as.numeric(format(bloom$Date, '%Y'))

## capture survey effort BEFORE filtering empty plant entries.
if (!'PlantGenusSpecies' %in% names(veg)) {
  veg$PlantGenusSpecies <- fix.white.space(paste(
    ifelse(is.na(veg$PlantGenus),'',veg$PlantGenus),
    ifelse(is.na(veg$PlantSpecies),'',veg$PlantSpecies),
    ifelse(is.na(veg$PlantVar),'',veg$PlantVar)))
}

##  In these standardized surveys, a quadrat containing no recorded
## flowers can be ABSENT from raw quadratveg.csv.  Therefore observed
## quadrat rows cannot define the floral denominator.  The full
## protocol design is reconstructed below (20 in 2012; 27 later).
## Preserve the four event/quadrat keys, even when that quadrat has no
## identified flowers. Do not infer survey effort from the
## plant-filtered rows below.
recorded.quads <- unique(veg[, c('Site', 'Year', 'SampleRound', 'Quadrat'),
                             drop=FALSE])
veg <- veg[!is.na(veg$PlantGenusSpecies) &
             nzchar(veg$PlantGenusSpecies),,drop=FALSE]

if (!'PlantGenusSpecies' %in% names(bloom)) {
  bloom$PlantGenusSpecies <- fix.white.space(paste(
    ifelse(is.na(bloom$PlantGenus),'',bloom$PlantGenus),
    ifelse(is.na(bloom$PlantSpecies),'',bloom$PlantSpecies),
    ifelse(is.na(bloom$PlantVar),'',bloom$PlantVar)))
}
bloom <- bloom[bloom$PlantGenusSpecies != "",]

## All veg and bloom plant names were reconciled against the same
## GBIF mapping as insect visit plants in relational/1prep.R. Do not
## normalize from a separate CSV here: that can discard morphospecies.

##  Protocol quadrats, rather than the subset
## with floral rows, define quadrat-level zeros and per-site means.
## 2012: 20 standard labels, later years: 27 standard labels.
quads.2012.design <- c(paste0(rep(c('C','E','N','S','W','ME'), each=3), 1:3),
                       'MN1', 'MN2')
quads.later.design <- paste0(
  rep(c('C','E','N','S','W','ME','MN','MS','MW'), each=3), 1:3)
## Excluded unknown labels are exported for manual correction, not silently
## admitted as a 28th or 29th quadrat. Their observations are not used in
## the standardized floral mean until the label is resolved.
valid.quadrat <- function(y, q) {
  ifelse(y == 2012, q %in% quads.2012.design,
         q %in% quads.later.design)
}
invalid.quad <- !valid.quadrat(veg$Year, as.character(veg$Quadrat)) |
                is.na(veg$Quadrat)
invalid.quad[is.na(invalid.quad)] <- TRUE
dir.create('../../skyIslands_saved/data/checks',
           recursive=TRUE, showWarnings=FALSE)
write.csv(veg[invalid.quad, ,drop=FALSE],
 file='../../skyIslands_saved/data/checks/nonstandard_vegetation_quadrats.csv',
 row.names=FALSE)
if (any(invalid.quad)) warning('## UPDATED [USER-QUADRATS]: ',
 sum(invalid.quad), ' floral rows have nonstandard quadrat IDs; ',
 'see checks/nonstandard_vegetation_quadrats.csv for review.')
veg <- veg[!invalid.quad, ,drop=FALSE]
## Survey rounds are identified from actual floral records and from
## documented Net effort: a fully empty floral round can have NO veg rows.
event.conditions <- read.csv(paste0(
 '../../skyIslands_saved/data/relational/relational/traditional/',
 'conditions-specimens-complete.csv'), stringsAsFactors=FALSE,
 colClasses=c(SampleRound='character'))
event.net <- event.conditions[!is.na(event.conditions$Method) &
  tolower(trimws(event.conditions$Method)) == 'net',
  c('Site','Year','SampleRound'), drop=FALSE]
event.net$Year <- as.numeric(as.character(event.net$Year))
recorded.rounds <- unique(recorded.quads[, c('Site','Year','SampleRound')])
recorded.rounds$SampleRound <- as.character(recorded.rounds$SampleRound)
recorded.rounds$Year <- as.numeric(as.character(recorded.rounds$Year))
event.net$SampleRound <- as.character(event.net$SampleRound)
## Do not turn vegetation without
## verified event geography into a newly inferred sampling location.
missing.veg.events <- dplyr::anti_join(recorded.rounds,
 event.net, by=c('Site','Year','SampleRound'))
write.csv(missing.veg.events,
 file='../../skyIslands_saved/data/checks/vegetation_without_verified_event.csv',
 row.names=FALSE)
if (nrow(missing.veg.events)) warning('## UPDATED: ',
 nrow(missing.veg.events),
 ' vegetation rounds lack verified site/round conditions; see checks.')
rounds <- unique(event.net)
rounds <- rounds[stats::complete.cases(rounds),,drop=FALSE]
quadrat.design <- data.frame(
 Year=sort(unique(rounds$Year)), stringsAsFactors=FALSE)
quadrat.design <- do.call(rbind, lapply(quadrat.design$Year, function(y) {
 data.frame(Year=y, Quadrat=if (y == 2012) quads.2012.design
            else quads.later.design, stringsAsFactors=FALSE)
}))
combos <- dplyr::inner_join(rounds, quadrat.design,
                             by='Year', relationship='many-to-many')
recorded.quads$Year <- as.numeric(as.character(recorded.quads$Year))
recorded.quads$SampleRound <- as.character(recorded.quads$SampleRound)
recorded.quads$RecordedVegetation <- TRUE
combos <- dplyr::left_join(combos,
  recorded.quads[, c('Site','Year','SampleRound','Quadrat',
                    'RecordedVegetation')],
  by=c('Site','Year','SampleRound','Quadrat'),
  relationship='one-to-one')
combos$RecordedVegetation[is.na(combos$RecordedVegetation)] <- FALSE
combos$ImputedZeroQuadrat <- !combos$RecordedVegetation
check_unique_key(combos, c('Site', 'Year', 'SampleRound', 'Quadrat'),
                 'Protocol vegetation quadrats')

## blooming flowers
veg.blooming <- veg[!is.na(veg$NumBlooms) &
                    veg$NumBlooms != '' & veg$NumBlooms != '0',]

## didn't count the exact number of flowers in 2017-2018, use
## midpoints of bins
veg.blooming$NumBloomsNum <- veg.blooming$NumBlooms
##  Use arithmetic midpoints for each bin.
## '<10' uses 5 (for the positive-count interval 1-9); others use 55,
## 550, not the former approximations 50 and 100.
veg.blooming$NumBloomsNum[veg.blooming$NumBlooms == "<10"] <- 5
veg.blooming$NumBloomsNum[veg.blooming$NumBlooms == "10-100"] <- 55
veg.blooming$NumBloomsNum[veg.blooming$NumBlooms == "100-1000"] <- 550

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

combos <- dplyr::left_join(combos, veg.bloom.sum.quad,
  by=c('Site', 'Year', 'SampleRound', 'Quadrat'),
  relationship='one-to-one')

cols.to.fill <-  c("FloralAbundance", "FloralRichness",
                   "FloralDiversity", "FloweringPlantAbundance",
                   "FloweringPlantDiversity")
## Missing protocol-quadrat records mean zero flowers under the stated
## sampling design, and enter the mean.
combos[, cols.to.fill][is.na(combos[, cols.to.fill])] <- 0
write.csv(combos, file='../data/veg_quadrat_effort.csv', row.names=FALSE)

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

#write.csv(veg.bloom.quad.sp, file="../data/veg.bloom.quad.sp.csv", row.names=FALSE)
## Site-level summaries are EVENT anchored: empty Net rounds appear
## with geographic/condition/floral data, and NA for specimen-derived
## site statistics (not imaginary specimen rows).
site.sum$Year <- as.numeric(as.character(site.sum$Year))
site.sum$SampleRound <- as.character(site.sum$SampleRound)
veg.bloom.mean$SampleRound <- as.character(veg.bloom.mean$SampleRound)
check_unique_key(veg.bloom.mean, keys.sr, 'Vegetation sample-round stats')
site.sum <- dplyr::left_join(rounds, site.sum, by=keys.sr,
                             relationship='one-to-one')
site.sum <- dplyr::left_join(site.sum, veg.bloom.mean, by=keys.sr,
                             relationship='one-to-one')
## Geographic and weather attributes are derived from distinct
## condition rows, not specimens, preventing capture-weighted site
## coordinates.
event.geo <- event.conditions %>%
  dplyr::filter(tolower(trimws(Method)) == 'net') %>%
  dplyr::mutate(Year=as.numeric(as.character(Year)),
                SampleRound=as.character(SampleRound)) %>%
  dplyr::distinct(ConditionsPK, .keep_all=TRUE) %>%
  dplyr::group_by(Site, Year, SampleRound) %>%
  dplyr::summarise(
    GeoLat=mean.unique(Lat), GeoLong=mean.unique(Long),
    GeoElev=mean.unique(Elev), GeoArea=mean.unique(Area),
    ## Carry environmental conditions into otherwise
    ## empty rows so a no-capture event is still fully described.
    EventMeanTempStart=mean.unique(TempStart),
    EventMeanTempEnd=mean.unique(TempEnd),
    EventMeanWindStart=mean.unique(WindStart),
    EventMeanWindEnd=mean.unique(WindEnd),
    EventAPi=mean.unique(APi), EventGDD=mean.unique(GDD),
    EventMeanTmean=mean.unique(MeanTmean),
    EventMeanVPD=mean.unique(MeanVPD),
    .groups='drop')
site.sum <- dplyr::left_join(site.sum, event.geo,
  by=keys.sr, relationship='one-to-one')
## The site-level relational geography is available from conditions
## even when an event has no specimens. Keep the familiar
## Lat/Elev/Area fields populated in sitestats.csv for those blank
## events.  Do not infer specimen abundance or richness from these
## event rows.
site.sum$Lat <- dplyr::coalesce(site.sum$Lat, site.sum$GeoLat)
site.sum$Elev <- dplyr::coalesce(site.sum$Elev, site.sum$GeoElev)
site.sum$Area <- dplyr::coalesce(site.sum$Area, site.sum$GeoArea)
write.csv(site.sum, file='../data/sitestats.csv', row.names=FALSE)

## floral richness across the entire meadow across sampling rounds
veg.year.sum <- veg.blooming %>%
    group_by(Site, Year) %>%
    summarise(TotalFloralRichness= length(unique(PlantGenusSpecies)))

write.csv(veg.year.sum, file="../data/veg_species_richness.csv",
          row.names=FALSE)

## *******************************************************************
## merging site and spec data
## *******************************************************************

## Merging specimen data with site characteristic data.
print("Before merge with site characteristics")
print(dim(spec.net))
## full join adds synthetic rows that are NOT specimens.
spec.net$Year <- as.numeric(as.character(spec.net$Year))
spec.net$SampleRound <- as.character(spec.net$SampleRound)
spec.net.before <- spec.net
spec.net <- dplyr::left_join(spec.net, site.sum, by=keys.sr,
                             relationship='many-to-one')
check_rows_preserved(spec.net.before, spec.net, 'Site stats -> specimens')
print("After merge with site characteristics")
print(dim(spec.net))

## Load the csv with the parasitism rate by site/SR/year/spp. Create a
## key that matches the spec.net dataset and then merge the species
## level information to the full dataset.
spec.indiv <- read_csv('../data/spstats_net.csv', show_col_types=FALSE)
## use the four natural key columns rather than pasted IDs.
network.keys <- c('Site', 'Year', 'SampleRound', 'GenusSpecies')
spec.indiv$Year <- as.numeric(as.character(spec.indiv$Year))
spec.indiv$SampleRound <- as.character(spec.indiv$SampleRound)

## The species summary includes unidentified netted specimens (blank
## GenusSpecies). read_csv() interprets blank strings as NA.  Preserve
## those specimen records, but do not treat 'unidentified' as a
## species when attaching species-level parasite and abundance
## summaries.  First check the actual sample-event keys; missing
## geography/round is a different problem and must not be silently
## dropped.
missing.event.key <- is.na(spec.indiv$Site) |
  !nzchar(trimws(as.character(spec.indiv$Site))) |
  is.na(spec.indiv$Year) |
  is.na(spec.indiv$SampleRound) |
  !nzchar(trimws(as.character(spec.indiv$SampleRound)))
if (any(missing.event.key)) {
  dir.create('../../skyIslands_saved/data/checks', recursive=TRUE,
             showWarnings=FALSE)
  write.csv(spec.indiv[missing.event.key, ,drop=FALSE],
            '../../skyIslands_saved/data/checks/species_summary_missing_event_keys.csv',
            row.names=FALSE)
  stop('Species-level summary has ', sum(missing.event.key),
       ' rows without Site/Year/SampleRound. See saved/data/checks.')
}

## Exclude unidentified *summary rows from this join
## only*, not the underlying specimens or their site-level abundances.
## Keep an list of those groups in skyIslands_saved.
missing.species <- is.na(spec.indiv$GenusSpecies) |
  !nzchar(trimws(as.character(spec.indiv$GenusSpecies)))
if (any(missing.species)) {
  dir.create('../../skyIslands_saved/data/checks', recursive=TRUE,
             showWarnings=FALSE)
  write.csv(spec.indiv[missing.species, ,drop=FALSE],
            '../../skyIslands_saved/data/checks/unidentified_net_species_summaries.csv',
            row.names=FALSE)
  message('## UPDATED: ', sum(missing.species),
          ' unidentified species-summary groups retained, ',
          'not matched as taxonomic species.')
}
spec.indiv.join <- spec.indiv[!missing.species, ,drop=FALSE]
check_unique_key(spec.indiv.join, network.keys, 'Identified species-level summary')

## Merging species data with individual-level data.
print('Before merge with ind level data')
print(dim(spec.net))

spec.net.before <- spec.net
## Join only identified species summaries; all specimen rows are
## retained by the left join and checked below.
spec.net <- dplyr::left_join(spec.net, spec.indiv.join, by=network.keys,
                             relationship='many-to-one')
check_rows_preserved(spec.net.before, spec.net, 'Species stats -> specimens')
print("After merge with ind level data")
print(dim(spec.net))

spec.net$GenusSpecies[is.na(spec.net$GenusSpecies)] <- ""
save(spec.net, file="../data/spec_net.Rdata")

