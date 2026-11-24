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

spec.all.methods <-
    read.csv('data/relational/relational/traditional/specimens-complete.csv',
             stringsAsFactors=FALSE, colClasses=c(SampleRound='character'))

setwd("../skyIslands/dataPrep")
source("src/misc.R")
source("src/prepParasites.R")
source("src/join_checks.R")
source("src/prepNets.R")
source("src/guild_networks.R")
source("src/specialization.R")
source("src/SpCalcMetrics.R")
source("src/vaznull2.R")

spec.all.methods$GenusSpecies <- fix.white.space(paste(
  ifelse(is.na(spec.all.methods$Genus), '', spec.all.methods$Genus),
  ifelse(is.na(spec.all.methods$Species), '', spec.all.methods$Species),
  ifelse(is.na(spec.all.methods$SubSpecies), '', spec.all.methods$SubSpecies)))

if (!'PlantGenusSpecies' %in% names(spec.all.methods)) {
  spec.all.methods$PlantGenusSpecies <- fix.white.space(paste(
    ifelse(is.na(spec.all.methods$PlantGenus), '', spec.all.methods$PlantGenus),
    ifelse(is.na(spec.all.methods$PlantSpecies), '', spec.all.methods$PlantSpecies),
    ifelse(is.na(spec.all.methods$PlantVar), '', spec.all.methods$PlantVar)))
}

missing.plant <- is.na(spec.all.methods$PlantGenusSpecies) |
                 spec.all.methods$PlantGenusSpecies == ''
if ('SourcePlantName' %in% names(spec.all.methods)) {
  spec.all.methods$PlantGenusSpecies[missing.plant] <-
    spec.all.methods$SourcePlantName[missing.plant]
}

spec.all.methods$Date <- as.Date(spec.all.methods$Date, format='%Y-%m-%d')
spec.all.methods$Doy <- as.numeric(strftime(spec.all.methods$Date, format='%j'))
spec.all.methods$Year <- as.numeric(format(spec.all.methods$Date,'%Y'))

## All plant names already come from the reconciled relational plant
## table; no second taxonomy pass is applied during specimen prep.

spec.all.methods$Int <- fix.white.space(paste(spec.all.methods$GenusSpecies,
                                  spec.all.methods$PlantGenusSpecies))
check_unique_key(spec.all.methods, 'UniqueID', 'Specimens loaded for analysis')

## *****************************************************************
## specimen-level parasite calculations
## *****************************************************************
## General output setup (also used by network and site-level analyses).
dir.create('../data/', showWarnings=FALSE)
dir.create('../data/networks', showWarnings=FALSE)
dir.create('../data/splevel_network_metrics', showWarnings=FALSE)

## Process assays and export the no-ID QC file using src/prepParasites.R.
parasite.prep <- prepareParasites(spec.all.methods)
spec.all.methods <- parasite.prep$spec
parasites <- parasite.prep$parasites  # Used by the parasite-network section.
rm(parasite.prep)

## Archive all methods after assay preparation, without method-specific
## summaries. This file is for reference, not used in downstream analyses.
write.csv(spec.all.methods, file='../data/spec_all_methods.csv', row.names=FALSE)

## All-method specimens were archived above. Analyses below use only Net.
spec.net <- spec.all.methods[
  !is.na(spec.all.methods$Method) & spec.all.methods$Method == 'Net',
  , drop=FALSE]
rm(spec.all.methods)  # Avoid accidentally mixing all-method and Net specimens.

## *****************************************************************
## Net sampling events and environmental conditions
## *****************************************************************

site.sr.keys <- c('Site', 'Year', 'SampleRound')
spec.net$SampleRound <- as.character(spec.net$SampleRound)

event.conditions <- read.csv(
  '../../skyIslands_saved/data/relational/relational/traditional/conditions-specimens-complete.csv',
  stringsAsFactors=FALSE, colClasses=c(SampleRound='character'))

## Ignore missing/nonfinite values; return NA for entirely missing groups.
mean.observed <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}
mean.unique <- function(x) mean.observed(unique(x))

## One record per condition, independent of the number of captures.
net.events <- event.conditions %>%
  filter(tolower(trimws(Method)) == 'net') %>%
  mutate(Year=as.numeric(as.character(Year)),
         SampleRound=as.character(SampleRound),
         Date=as.Date(Date, format='%Y-%m-%d'),
         Doy=as.numeric(format(Date, '%j')))

## A missing condition ID or sampling key cannot define a verified event.
check.event.cols <- c('ConditionsPK', site.sr.keys)
invalid.event <- !complete.cases(net.events[, check.event.cols]) |
  Reduce(`|`, lapply(net.events[, check.event.cols], function(x)
    !nzchar(trimws(as.character(x)))))
if (any(invalid.event)) {
  stop('Net conditions contain missing IDs or sampling keys; inspect event.conditions.')
}
if (!nrow(net.events)) stop('No documented Net sampling events.')
net.events <- net.events %>% distinct(ConditionsPK, .keep_all=TRUE)

## One authoritative row per site-round, including zero-capture rounds.
## Weather and day of year are averaged across distinct condition records.
net.site.sr.events <- net.events %>%
  group_by(Site, Year, SampleRound) %>%
  summarise(
    Lat=mean.unique(Lat), Long=mean.unique(Long),
    Elev=mean.unique(Elev), Area=mean.unique(Area),
    SRDoy=mean.observed(Doy),
    EventMeanTempStart=mean.observed(TempStart),
    EventMeanTempEnd=mean.observed(TempEnd),
    EventMeanWindStart=mean.observed(WindStart),
    EventMeanWindEnd=mean.observed(WindEnd),
    EventAPi=mean.unique(APi), EventGDD=mean.unique(GDD),
    EventMeanTmean=mean.unique(MeanTmean),
    EventMeanVPD=mean.unique(MeanVPD),
    across(c(Tmean_anom, Tmax_anom, Tmin_anom, Precip_anom, VPD_anom),
           mean.unique),
    .groups='drop')
check_unique_key(net.site.sr.events, site.sr.keys, 'Net sampling events')

## Do not silently lose observed specimen rounds when anchoring to events.
net.rounds.without.events <- spec.net %>%
  select(all_of(site.sr.keys)) %>% distinct() %>%
  anti_join(net.site.sr.events, by=site.sr.keys)
if (nrow(net.rounds.without.events)) {
  warning(nrow(net.rounds.without.events),
          ' specimen rounds lack verified Net conditions; inspect net.rounds.without.events.')
}

## ***********************************************************************
## Net specimen summaries at Site x Year x SampleRound and species grains
## Naming: site = Site, year = Year, sr = SampleRound, sp = species,
## quad = Quadrat; 'sum' indicates grouped data, not individual records.
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


calcNetSummaryStats <- function(net.specimens){
    net.site.sr.sp.sum <- net.specimens %>%
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
    net.site.sr.sum <- net.specimens %>%
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
                  CrithidiaParasitismRate=mean(CrithidiaPresence, na.rm=TRUE),
                  ApicystisParasitismRate=mean(ApicystisSpp, na.rm=TRUE),
                  CrithidiaBombusParasitismRate= mean(CrithidiaPresence[Genus == 'Bombus'],
                                                      na.rm=TRUE),
                  CrithidiaHBParasitismRate= mean (CrithidiaPresence
                                                   [GenusSpecies == 'Apis mellifera'],
                                                   na.rm=TRUE))

    net.site.year.sp.sum <- net.specimens %>%
        group_by(Site, Year, GenusSpecies, Genus) %>%
        summarise(Abundance = length(GenusSpecies))


    net.site.sr.sp.counts <- net.specimens %>%
      group_by(Site, Year, SampleRound, GenusSpecies, Genus) %>%

      summarise(AbundanceSYR = length(GenusSpecies))

    net.site.year.bombus.sum <- net.site.year.sp.sum[
      net.site.year.sp.sum$Genus == 'Bombus', , drop=FALSE]
    net.site.year.bombus.sum$Genus  <- NULL

    ## Export the same legacy filenames so downstream scripts keep working.
    write.csv(net.site.year.bombus.sum,
              file='../data/bombus_year_site_net.csv', row.names=FALSE)
    write.csv(net.site.year.sp.sum,
              file='../data/sp_year_site_net.csv', row.names=FALSE)
    write.csv(net.site.sr.sp.sum,
              file='../data/spstats_net.csv', row.names=FALSE)
    write.csv(net.site.sr.sp.counts,
              file='../data/sp_year_site_round_net.csv', row.names=FALSE)
    write.csv(net.site.sr.sp.counts,
              file='../data/sp_year_site_round.csv', row.names=FALSE)
    return(net.site.sr.sum)
}

net.site.sr.sum <- calcNetSummaryStats(spec.net)

## ***********************************************************************
## NETTED genus-by-year and ID-completeness diagnostics.
## Export diagnostic CSVs and plots to skyIslands_saved/data/checks.
## Source is Net specimens before joins, so each row is counted once.
## ***********************************************************************
dir.guild.audit <- '../../skyIslands_saved/data/checks/guild_id_review'
dir.create(dir.guild.audit, recursive=TRUE, showWarnings=FALSE)
source("src/guild_review.R")

table(net.site.sr.sum$Site)

## Site x Year x SampleRound: Net capture statistics only.
check_unique_key(net.site.sr.sum, site.sr.keys, 'Net sample-round stats')


## ***********************************************************************
## Parasite summaries at species, genus and site-round grains
## ***********************************************************************

parasite.site.sr.sp.sum <- spec.net %>%
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

parasite.site.sr.gen.sum <- spec.net %>%
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

parasite.site.sr.sum <- spec.net %>%
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



## Keep the saved object name `all.sums` for existing RData consumers.
all.sums <- left_join(parasite.site.sr.sp.sum, parasite.site.sr.gen.sum,
  by=c('Site', 'Year', 'SampleRound', 'Genus'),
  relationship='many-to-one')
all.sums <- left_join(all.sums, parasite.site.sr.sum,
  by=c('Site', 'Year', 'SampleRound'),
  relationship='many-to-one')
all.sums$Project <- "SI"

save(all.sums, file="../data/SI_parasite_sums.Rdata")


## *******************************************************************
## Dataset-wide plant-pollinator specialization (Net specimens)
## *******************************************************************
## Exclude empty plant or pollinator IDs from observed interactions.
## Apply the same basic ID consistency checks before calculating
## dataset-wide plant-pollinator traits and event-level guild networks.
## Raw specimen records themselves are never deleted.
net.specimens.for.interactions <- audit_network_ids(spec.net)$usable

net.interactions.by.species.pair <- aggregate(list(abund=net.specimens.for.interactions$GenusSpecies),
                      list(GenusSpecies=net.specimens.for.interactions$GenusSpecies,
                           PlantGenusSpecies=net.specimens.for.interactions$PlantGenusSpecies),
                      length)

net.plant.pollinator.matrix <- samp2site.spp(
  net.interactions.by.species.pair$PlantGenusSpecies,
  net.interactions.by.species.pair$GenusSpecies,
  net.interactions.by.species.pair$abund, FUN=sum)

net.plant.pollinator.network.traits <- specieslevel(net.plant.pollinator.matrix)
## Rarefied degree for plants and pollinators
rare.plants.degree <- apply(net.plant.pollinator.matrix, 1, chao1)
rare.pols.degree <- apply(net.plant.pollinator.matrix, 2, chao1)

net.plant.pollinator.traits.table <- data.frame(
  GenusSpecies=unlist(sapply(net.plant.pollinator.network.traits, rownames)),
  do.call(rbind, net.plant.pollinator.network.traits))

net.plant.pollinator.traits.table$r.degree <- rare.pols.degree[match(
  net.plant.pollinator.traits.table$GenusSpecies, names(rare.pols.degree))]
missing.degree <- is.na(net.plant.pollinator.traits.table$r.degree)
net.plant.pollinator.traits.table$r.degree[missing.degree] <-
  rare.plants.degree[match(
    net.plant.pollinator.traits.table$GenusSpecies[missing.degree],
    names(rare.plants.degree))]

rownames(net.plant.pollinator.traits.table) <- NULL

write.csv(net.plant.pollinator.traits.table, file='../data/networks_traits.csv', row.names=FALSE)

## *******************************************************************
## Dataset-wide parasite-pollinator specialization (Net specimens)
## specialization etc. across all SI
## *******************************************************************
net.screened.specimens <- spec.net[spec.net$Apidae == 1 & !is.na(spec.net$Apidae),]

net.parasite.by.species <- aggregate(net.screened.specimens[, parasites],
                      list(GenusSpecies=net.screened.specimens$GenusSpecies),
                                              sum, na.rm=TRUE)

net.screened.by.species <- table(net.screened.specimens$GenusSpecies)
net.screened.by.species

net.screened.by.genus <- table(net.screened.specimens$Genus)
net.screened.by.genus

## proportion of individuals screened
net.parasite.by.species[, parasites] <- net.parasite.by.species[, parasites]/
    net.screened.by.species[net.parasite.by.species$GenusSpecies]

net.parasite.matrix <- net.parasite.by.species
rownames(net.parasite.matrix) <- net.parasite.matrix$GenusSpecies
net.parasite.matrix$GenusSpecies <- NULL
net.parasite.matrix <- as.matrix(net.parasite.matrix)


net.parasite.network.traits <- specieslevel(net.parasite.matrix)
## Parasite specialization statistics for screened Net specimens.

net.parasite.traits.table <- data.frame(
  GenusSpecies=unlist(sapply(net.parasite.network.traits, rownames)),
  do.call(rbind, net.parasite.network.traits))

rownames(net.parasite.traits.table) <- NULL

write.csv(net.parasite.traits.table, file='../data/parasitetraits.csv', row.names=FALSE)

## *******************************************************************
## Site x Year x SampleRound plant-pollinator and parasite networks
## *******************************************************************

## Construct an explicitly filtered named list of plant x insect matrices
## for every pollinator guild, at Site x Year x SampleRound only.
## The event export also records true netting rounds with zero observations.
## Reload the network helper when this block is run interactively after
## editing the file; otherwise a previous definition remains in memory.
source('src/guild_networks.R')
net.site.sr.guild.networks <- make_guild_networks(
  spec.net, net.site.sr.events %>% select(all_of(site.sr.keys)),
  guilds=pollinator_guilds(bee.families, butterfly.hawk.families),
  network_dir='../data/networks',
  metrics_dir='../data/splevel_network_metrics',
  audit_dir='../../skyIslands_saved/data/checks/network_id_review')
net.site.sr.bee.network.metrics <- net.site.sr.guild.networks$bee_metrics

## Parasites
net.parasite.assays.wide <- net.screened.specimens %>%
  dplyr::select(UniqueID, GenusSpecies, Site, Year, SampleRound, AscosphaeraSpp,
         ApicystisSpp, CrithidiaExpoeki, CrithidiaBombi,
         CrithidiaMellificae,
         NosemaBombi, NosemaCeranae)

net.parasite.assays.long <- net.parasite.assays.wide %>%
  pivot_longer(cols=c("AscosphaeraSpp",
                      "ApicystisSpp", "CrithidiaExpoeki",
                      "CrithidiaBombi", "CrithidiaMellificae",
                      "NosemaBombi", "NosemaCeranae"),
               names_to = "Parasite", values_to = "count")
net.parasite.assays.long <- as.data.frame(net.parasite.assays.long)

## an NA assay is unobserved, not an absent parasite edge;
## omit those particular assays before aggregating species prevalence.
net.parasite.assays.long <- net.parasite.assays.long[!is.na(net.parasite.assays.long$count), ]

net.site.sr.parasite.networks <- makeNets(net.parasite.assays.long, net.type="YearSR",
                           species=c("Pollinator",
                                     "Parasite"),
                           lower.level="GenusSpecies",
                           higher.level="Parasite",
                           poll.groups="Bees")


## Annual mean parasite networks are intentionally not generated.
## Each sampling round's network is retained in net.site.sr.parasite.networks.

## Attach site-round plant-pollinator network metrics to Net specimens
print("specimen dim before network metric year-sr merge")
dim(spec.net)

net.site.sr.bee.network.metrics <- net.site.sr.bee.network.metrics[
  net.site.sr.bee.network.metrics$speciesType == 'higher.level', , drop=FALSE]
net.site.sr.bee.network.metrics$SpSiteYear <- NULL
site.sr.sp.keys <- c('Site', 'Year', 'SampleRound', 'GenusSpecies')
net.site.sr.bee.network.metrics$Year <- as.numeric(as.character(net.site.sr.bee.network.metrics$Year))
net.site.sr.bee.network.metrics$SampleRound <- as.character(net.site.sr.bee.network.metrics$SampleRound)
spec.net$SampleRound <- as.character(spec.net$SampleRound)
check_unique_key(net.site.sr.bee.network.metrics, site.sr.sp.keys, 'Bee network metrics')
spec.net.before.join <- spec.net
spec.net <- dplyr::left_join(spec.net, net.site.sr.bee.network.metrics, by=site.sr.sp.keys,
                             relationship='many-to-one')
check_rows_preserved(spec.net.before.join, spec.net, 'Network metrics -> specimens')
dim(spec.net)
print("specimen dim after network metric year-sr merge")

spec.net$speciesType <- NULL

## *******************************************************************
## Net specimen data checks
## *******************************************************************

print(paste("Pollinator species", length(unique(spec.net$GenusSpecies))))
print(paste("Bee species", length(unique(spec.net$GenusSpecies[
                                                spec.net$Family %in%
                                                bee.families]))))
print(paste("Syrphid species", length(unique(spec.net$GenusSpecies[
                                                spec.net$Family == "Syrphidae"]))))
print(paste("Plant species", length(unique(spec.net$PlantGenusSpecies))))
print(paste("Pollinator genera", length(unique(spec.net$Genus))))
print(paste("Interactions", length(unique(spec.net$Int))))
print(paste("Specimens", nrow(spec.net)))

## *******************************************************************
## Vegetation and bloom preparation
## *******************************************************************

bloom <-
    read.csv("../../skyIslands_saved/data/relational/relational/traditional/bloom-complete.csv",
             stringsAsFactors=FALSE)
veg <-
    read.csv("../../skyIslands_saved/data/relational/relational/traditional/veg-complete.csv",
             stringsAsFactors=FALSE, colClasses=c(SampleRound='character'))

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
veg.site.sr.quads.recorded <- unique(veg[, c('Site', 'Year', 'SampleRound', 'Quadrat'),
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
## Verify vegetation rounds against documented Net events. A zero-bloom
## round may have no vegetation observations; Net effort anchors the design.
veg.site.sr.rounds.recorded <- unique(veg.site.sr.quads.recorded[
  , c('Site', 'Year', 'SampleRound')])
veg.site.sr.rounds.recorded$SampleRound <- as.character(veg.site.sr.rounds.recorded$SampleRound)
veg.site.sr.rounds.recorded$Year <- as.numeric(as.character(veg.site.sr.rounds.recorded$Year))
## Do not turn vegetation without
## verified event geography into a newly inferred sampling location.
veg.rounds.without.net.events <- dplyr::anti_join(veg.site.sr.rounds.recorded,
 net.site.sr.events, by=site.sr.keys)
write.csv(veg.rounds.without.net.events,
 file='../../skyIslands_saved/data/checks/vegetation_without_verified_event.csv',
 row.names=FALSE)
if (nrow(veg.rounds.without.net.events)) warning('## UPDATED: ',
 nrow(veg.rounds.without.net.events),
 ' vegetation rounds lack verified site/round conditions; see checks.')
veg.quad.labels.by.year <- data.frame(
 Year=sort(unique(net.site.sr.events$Year)), stringsAsFactors=FALSE)
veg.quad.labels.by.year <- do.call(rbind, lapply(veg.quad.labels.by.year$Year, function(y) {
 data.frame(Year=y, Quadrat=if (y == 2012) quads.2012.design
            else quads.later.design, stringsAsFactors=FALSE)
}))
veg.site.sr.quads.design <- net.site.sr.events %>%
  select(all_of(site.sr.keys)) %>%
  inner_join(veg.quad.labels.by.year, by='Year', relationship='many-to-many')
veg.site.sr.quads.recorded$Year <- as.numeric(as.character(veg.site.sr.quads.recorded$Year))
veg.site.sr.quads.recorded$SampleRound <- as.character(veg.site.sr.quads.recorded$SampleRound)
veg.site.sr.quads.recorded$RecordedVegetation <- TRUE
veg.site.sr.quads.design <- dplyr::left_join(veg.site.sr.quads.design,
  veg.site.sr.quads.recorded[, c('Site','Year','SampleRound','Quadrat',
                    'RecordedVegetation')],
  by=c('Site','Year','SampleRound','Quadrat'),
  relationship='one-to-one')
veg.site.sr.quads.design$RecordedVegetation[is.na(veg.site.sr.quads.design$RecordedVegetation)] <- FALSE
veg.site.sr.quads.design$ImputedZeroQuadrat <- !veg.site.sr.quads.design$RecordedVegetation
check_unique_key(veg.site.sr.quads.design, c('Site', 'Year', 'SampleRound', 'Quadrat'),
                 'Protocol vegetation quadrats')

## Filter to positive blooming-plant observations.
veg.blooming <- veg[!is.na(veg$NumBlooms) &
                    veg$NumBlooms != '' & veg$NumBlooms != '0',]

## Use bin midpoints for years 2017-2018, which recorded ranges.
veg.blooming$NumBloomsNum <- veg.blooming$NumBlooms
##  Use arithmetic midpoints for each bin.
## '<10' uses 5 (for the positive-count interval 1-9); others use 55,
## 550, not the former approximations 50 and 100.
veg.blooming$NumBloomsNum[veg.blooming$NumBlooms == "<10"] <- 5
veg.blooming$NumBloomsNum[veg.blooming$NumBlooms == "10-100"] <- 55
veg.blooming$NumBloomsNum[veg.blooming$NumBlooms == "100-1000"] <- 550

veg.blooming$NumBloomsNum <- as.numeric(veg.blooming$NumBloomsNum)

## *******************************************************************
## Vegetation summaries from quadrats
## *******************************************************************

## Species-by-quadrat floral abundance precedes the site-round means.
veg.site.sr.quad.sp.sum <- veg.blooming %>%
    group_by(Quadrat, Site, Year, PlantGenusSpecies, SampleRound) %>%
    summarise(FloralAbundance = sum(NumBloomsNum),
              FloweringPlantAbundance= sum(PlantCount, na.rm=TRUE))

veg.site.sr.quad.sp.sum$KeyQuad <- paste(veg.site.sr.quad.sp.sum$Site,
                              veg.site.sr.quad.sp.sum$Quadrat,
                              veg.site.sr.quad.sp.sum$Year,
                              veg.site.sr.quad.sp.sum$SampleRound)

veg.quad.floral.shannon <- tapply(veg.site.sr.quad.sp.sum$FloralAbundance,
                  veg.site.sr.quad.sp.sum$KeyQuad,
                  vegan::diversity)

veg.quad.plant.shannon <- tapply(veg.site.sr.quad.sp.sum$FloweringPlantAbundance,
                  veg.site.sr.quad.sp.sum$KeyQuad,
                  vegan::diversity)


veg.site.sr.quad.sum <- veg.site.sr.quad.sp.sum %>%
    group_by(Quadrat, Site, Year, SampleRound) %>%
    summarise(FloralAbundance = sum(FloralAbundance),
              FloralRichness= length(unique(PlantGenusSpecies)),
              FloweringPlantAbundance = sum(FloweringPlantAbundance)
              )

veg.quad.lookup.keys <-   names(veg.quad.floral.shannon)
names(veg.quad.floral.shannon) <- NULL

## Attach quadrat Shannon indices to quadrat summaries.
veg.site.sr.quad.sum$FloralDiversity <-  veg.quad.floral.shannon[match(
    paste(veg.site.sr.quad.sum$Site,
          veg.site.sr.quad.sum$Quadrat,
          veg.site.sr.quad.sum$Year,
          veg.site.sr.quad.sum$SampleRound),
    veg.quad.lookup.keys)]

veg.site.sr.quad.sum$FloweringPlantDiversity <-
    veg.quad.plant.shannon[match(
        paste(veg.site.sr.quad.sum$Site,
              veg.site.sr.quad.sum$Quadrat,
              veg.site.sr.quad.sum$Year,
              veg.site.sr.quad.sum$SampleRound),
        veg.quad.lookup.keys)]

## Join quadrat-level floral statistics to the full protocol quadrat design

veg.site.sr.quads.design <- dplyr::left_join(veg.site.sr.quads.design, veg.site.sr.quad.sum,
  by=c('Site', 'Year', 'SampleRound', 'Quadrat'),
  relationship='one-to-one')

veg.quad.zero.columns <-  c("FloralAbundance", "FloralRichness",
                   "FloralDiversity", "FloweringPlantAbundance",
                   "FloweringPlantDiversity")
## Missing protocol-quadrat records mean zero flowers under the stated
## sampling design, and enter the mean.
veg.site.sr.quads.design[, veg.quad.zero.columns][
  is.na(veg.site.sr.quads.design[, veg.quad.zero.columns])] <- 0
write.csv(veg.site.sr.quads.design, file='../data/veg_quadrat_effort.csv', row.names=FALSE)

veg.site.sr.sum <- veg.site.sr.quads.design %>%
    group_by(Site, Year, SampleRound) %>%
    summarise(MeanFloralAbundance = mean(FloralAbundance),
              MeanFloralRichness= mean(FloralRichness),
              MeanFloralDiversity =mean(FloralDiversity),
              MeanFloweringPlantDiversity =mean(FloweringPlantDiversity),
              MeanFloweringPlantAbundance =mean(FloweringPlantAbundance)
              )

veg.site.sr.sum[order(veg.site.sr.sum$MeanFloralRichness), ]

write.csv(veg.site.sr.sum, file="../data/veg.csv", row.names=FALSE)

#write.csv(veg.site.sr.quad.sp.sum,
#          file='../data/veg.bloom.quad.sp.csv', row.names=FALSE)
## Combine verified Net events, captures and independent vegetation summaries.
check_unique_key(veg.site.sr.sum, site.sr.keys, 'Vegetation sample-round stats')
net.site.sr.sum <- net.site.sr.events %>%
  left_join(net.site.sr.sum, by=site.sr.keys, relationship='one-to-one') %>%
  left_join(veg.site.sr.sum, by=site.sr.keys, relationship='one-to-one')

## Only verified zero-capture rounds receive zeros for capture metrics.
## Parasite rates remain NA, and vegetation and conditions are unchanged.
empty.capture <- is.na(net.site.sr.sum$PollAbundance)
capture.cols <- grep(
  '^(Poll|Bee|Syrphid|HB|Bombus|NonBombusHB|ButterflyHawk|VisitedFloral).*(Abundance|Richness|Diversity)$',
  names(net.site.sr.sum), value=TRUE)
for (col in capture.cols) net.site.sr.sum[[col]][empty.capture] <- 0
message('Documented Net rounds without specimens: ', sum(empty.capture))

write.csv(net.site.sr.sum, file='../data/sitestats.csv', row.names=FALSE)

## floral richness across the entire meadow across sampling rounds
veg.year.sum <- veg.blooming %>%
    group_by(Site, Year) %>%
    summarise(TotalFloralRichness= length(unique(PlantGenusSpecies)))

write.csv(veg.year.sum, file="../data/veg_species_richness.csv",
          row.names=FALSE)

## *******************************************************************
## Merge site-round and species summaries onto Net specimens
## *******************************************************************

## Attach Site x Year x SampleRound summaries to Net specimen rows.
print("Before merge with site characteristics")
print(dim(spec.net))
## Use a left join: no synthetic specimen rows should be created.
spec.net$Year <- as.numeric(as.character(spec.net$Year))
spec.net$SampleRound <- as.character(spec.net$SampleRound)
spec.net.before.join <- spec.net
## Keep original specimen conditions under their original names; overlapping
## site-round summaries receive an explicit .site.sr suffix.
spec.net <- dplyr::left_join(spec.net, net.site.sr.sum, by=site.sr.keys,
                             suffix=c('', '.site.sr'),
                             relationship='many-to-one')
check_rows_preserved(spec.net.before.join, spec.net, 'Site stats -> specimens')
print("After merge with site characteristics")
print(dim(spec.net))

## Load Site x Year x SampleRound x Species statistics from the Net
## export and attach them only to identified Net specimens.
net.site.sr.sp.stats <- read_csv('../data/spstats_net.csv', show_col_types=FALSE,
                                        col_types=cols(SampleRound=col_character()))
## Use the four natural key columns rather than pasted IDs.
site.sr.sp.keys <- c('Site', 'Year', 'SampleRound', 'GenusSpecies')
net.site.sr.sp.stats$Year <- as.numeric(as.character(net.site.sr.sp.stats$Year))
net.site.sr.sp.stats$SampleRound <- as.character(net.site.sr.sp.stats$SampleRound)

## The species summary includes unidentified netted specimens (blank
## GenusSpecies). read_csv() interprets blank strings as NA. Preserve
## those specimen records, but do not treat 'unidentified' as a
## species when attaching species-level parasite and abundance
## summaries.  First check the actual sample-event keys; missing
## geography/round is a different problem and must not be silently
## dropped.
missing.event.key <- is.na(net.site.sr.sp.stats$Site) |
  !nzchar(trimws(as.character(net.site.sr.sp.stats$Site))) |
  is.na(net.site.sr.sp.stats$Year) |
  is.na(net.site.sr.sp.stats$SampleRound) |
  !nzchar(trimws(as.character(net.site.sr.sp.stats$SampleRound)))
if (any(missing.event.key)) {
  dir.create('../../skyIslands_saved/data/checks', recursive=TRUE,
             showWarnings=FALSE)
  write.csv(net.site.sr.sp.stats[missing.event.key, ,drop=FALSE],
            '../../skyIslands_saved/data/checks/species_summary_missing_event_keys.csv',
            row.names=FALSE)
  stop('Species-level summary has ', sum(missing.event.key),
       ' rows without Site/Year/SampleRound. See saved/data/checks.')
}

## Exclude unidentified *summary rows from this join
## only*, not the underlying specimens or their site-level abundances.
## Keep a list of those groups in skyIslands_saved.
missing.species <- is.na(net.site.sr.sp.stats$GenusSpecies) |
  !nzchar(trimws(as.character(net.site.sr.sp.stats$GenusSpecies)))
if (any(missing.species)) {
  dir.create('../../skyIslands_saved/data/checks', recursive=TRUE,
             showWarnings=FALSE)
  write.csv(net.site.sr.sp.stats[missing.species, ,drop=FALSE],
            '../../skyIslands_saved/data/checks/unidentified_net_species_summaries.csv',
            row.names=FALSE)
  message('## UPDATED: ', sum(missing.species),
          ' unidentified species-summary groups retained, ',
          'not matched as taxonomic species.')
}
net.site.sr.sp.stats.join <- net.site.sr.sp.stats[
  !missing.species, , drop=FALSE]
check_unique_key(net.site.sr.sp.stats.join, site.sr.sp.keys,
                 'Identified species-level summary')

## Attach species-round statistics to Net specimen records.
print('Before merge with ind level data')
print(dim(spec.net))

spec.net.before.join <- spec.net
## Join only identified species summaries; all specimen rows are
## retained by the left join and checked below.
spec.net <- dplyr::left_join(spec.net, net.site.sr.sp.stats.join, by=site.sr.sp.keys,
                             relationship='many-to-one')
check_rows_preserved(spec.net.before.join, spec.net, 'Species stats -> specimens')
print("After merge with ind level data")
print(dim(spec.net))

spec.net$GenusSpecies[is.na(spec.net$GenusSpecies)] <- ""
save(spec.net, file="../data/spec_net.Rdata")

