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
source("src/guild_networks.R")

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

## Process assays and export the no-ID QC file using src/prepParasites.R.
parasite.prep <- prepareParasites(spec.all.methods)
spec.all.methods <- parasite.prep$spec
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

stopifnot(nrow(net.events) > 0,
          all(complete.cases(net.events[, c("ConditionsPK", site.sr.keys)])))
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
## Other moth families (e.g., Erebidae, Noctuidae) are not included
## since they were not identified to species consistently between years.

butterfly.hawk.families <- c('Hesperiidae', 'Lycaenidae', 'Nymphalidae',
                            'Papilionidae', 'Pieridae', 'Riodinidae',
                            'Sphingidae')
focal.pollinator.families <- c(bee.families, 'Syrphidae',
                               butterfly.hawk.families)
## Keep the existing Bees+Syrphids+Lep columns unchanged. Add a separate
## combined guild that includes bee flies (Bombyliidae).
all.pollinator.families <- c(focal.pollinator.families, 'Bombyliidae')


## Ecological summaries: no individuals is zero richness, but no
## identified taxa among captured individuals is unknown richness.
## Observed taxonomic richness
observed.richness <- function(x) {
  richness <- sum(table(x[x != '']) > 0)
  if (length(x) > 0 && richness == 0) NA_integer_ else richness
}

## Shannon diversity is undefined for empty communities
shannon.counts <- function(x) {
  if (isTRUE(sum(x) == 0)) NA_real_ else vegan::diversity(x)
}

## Shannon diversity from individual taxonomic observations
observed.shannon <- function(x) {
  shannon.counts(table(x[x != '']))
}

sum.observed <- function(x) {
  if (all(is.na(x))) NA_real_ else sum(x, na.rm=TRUE)
}

calcNetSummaryStats <- function(net.specimens){
    net.site.sr.sp.sum <- net.specimens %>%
        group_by(Site, Year, SampleRound, GenusSpecies) %>%
        summarise(Abundance = length(GenusSpecies),
                  SpParasitismRate = mean.observed(ParasitePresence),
                  SpCrithidiaParasitismRate = mean.observed(CrithidiaPresence),
                  SpApicystisParasitismRate = mean.observed(ApicystisSpp),
                  SpCrithidiaBombusParasitismRate = mean.observed(CrithidiaPresence[Genus %in% "Bombus"]),
                  SpCrithidiaHBParasitismRate = mean.observed(CrithidiaPresence
                                                   [GenusSpecies %in% "Apis mellifera"]),
                  SpParasitism = sum.observed(ParasitePresence),
                  SpCrithidiaPresence= sum.observed(CrithidiaPresence),
                  SpApicystisSpp = sum.observed(ApicystisSpp),
                  SpNosemaBombi= sum.observed(NosemaBombi),
                  SpNosemaCeranae = sum.observed(NosemaCeranae),
                  SpAscosphaeraSpp= sum.observed(AscosphaeraSpp),
                  SpCrithidiaExpoeki = sum.observed(CrithidiaExpoeki),
                  SpCrithidiaMellificae = sum.observed(CrithidiaMellificae),
                  SpCrithidiaBombi = sum.observed(CrithidiaBombi),
                  SpCrithidiaSpp = sum.observed(CrithidiaSpp),
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
                      length(GenusSpecies[Family %in% 'Syrphidae']),
                  BeeFlyAbundance = sum(Family %in% 'Bombyliidae'),
                  HBAbundance = sum(GenusSpecies %in% 'Apis mellifera'),
                  BombusAbundance = sum(Genus %in% 'Bombus'),
                  NonBombusHBAbundance =
                      sum(!is.na(Genus) & nzchar(trimws(Genus)) & !Genus %in% c('Bombus', 'Apis') &
                          Family %in% bee.families),
                  ## New Lepidoptera and combined counts.
                  ButterflyHawkAbundance =
                      sum(Family %in% butterfly.hawk.families),
                  BeeSyrphLepAbundance =
                      sum(Family %in% focal.pollinator.families),
                  BeeSyrphBeeFlyLepAbundance =
                      sum(Family %in% all.pollinator.families),

                  PollRichness= observed.richness(GenusSpecies),
                  BeeRichness= observed.richness(GenusSpecies[Family %in%
                                                          bee.families]),
                  SyrphidRichness= observed.richness(
                    GenusSpecies[Family %in% 'Syrphidae']),
                  BeeFlyRichness = observed.richness(
                    GenusSpecies[Family %in% 'Bombyliidae']),
                  BombusRichness= observed.richness(
                    GenusSpecies[Genus %in% 'Bombus']),
                  ## Lepidoptera and
                  ## combined taxonomic-label richness. 
                  ButterflyHawkRichness = observed.richness(
                      GenusSpecies[Family %in% butterfly.hawk.families]),
                  BeeSyrphLepRichness = observed.richness(
                      GenusSpecies[Family %in% focal.pollinator.families]),
                  BeeSyrphBeeFlyLepRichness = observed.richness(
                      GenusSpecies[Family %in% all.pollinator.families]),

                  PollDiversity=observed.shannon(GenusSpecies),
                  BeeDiversity=observed.shannon(GenusSpecies[Family %in% bee.families]),
                  SyrphidDiversity=observed.shannon(GenusSpecies[Family %in% 'Syrphidae']),
                  BeeFlyDiversity = observed.shannon(
                    GenusSpecies[Family %in% 'Bombyliidae']),
                  BombusDiversity=observed.shannon(GenusSpecies[Genus %in% 'Bombus']),
                  ## Shannon for three added groups.
                  ButterflyHawkDiversity = observed.shannon(
                      GenusSpecies[Family %in% butterfly.hawk.families]),
                  BeeSyrphLepDiversity = observed.shannon(
                      GenusSpecies[Family %in% focal.pollinator.families]),
                  BeeSyrphBeeFlyLepDiversity = observed.shannon(
                      GenusSpecies[Family %in% all.pollinator.families]),

                  VisitedFloralRichness= observed.richness(PlantGenusSpecies),
                  VisitedFloralDiversity=observed.shannon(PlantGenusSpecies),

                  SiteParasitismRate=mean.observed(ParasitePresence),
                  MeanParasiteRichness=mean.observed(ParasiteRichness),
                  HBSiteParasitismRate=mean.observed(
                      ParasitePresence[GenusSpecies %in% 'Apis mellifera']),
                  BombusSiteParasitismRate=mean.observed(
                      ParasitePresence[Genus %in% 'Bombus']),
                  CrithidiaParasitismRate=mean.observed(CrithidiaPresence),
                  ApicystisParasitismRate=mean.observed(ApicystisSpp),
                  CrithidiaBombusParasitismRate= mean.observed(CrithidiaPresence[Genus %in% 'Bombus']),
                  CrithidiaHBParasitismRate= mean.observed(CrithidiaPresence
                                                   [GenusSpecies %in% 'Apis mellifera']))

    net.site.year.sp.sum <- net.specimens %>%
        group_by(Site, Year, GenusSpecies, Genus) %>%
        summarise(Abundance = length(GenusSpecies))


    net.site.sr.sp.counts <- net.specimens %>%
      group_by(Site, Year, SampleRound, GenusSpecies, Genus) %>%

      summarise(AbundanceSYR = length(GenusSpecies))

    net.site.year.bombus.sum <- net.site.year.sp.sum[
      net.site.year.sp.sum$Genus %in% 'Bombus', , drop=FALSE]
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
    return(list(site=net.site.sr.sum, species=net.site.sr.sp.sum))
}

net.summary <- calcNetSummaryStats(spec.net)
net.site.sr.sum <- net.summary$site
net.site.sr.sp.stats <- net.summary$species

## ***********************************************************************
## NETTED genus-by-year and ID-completeness diagnostics.
## Export diagnostic CSVs and plots to skyIslands_saved/data/checks.
## Source is Net specimens before joins, so each row is counted once.
## ***********************************************************************
dir.guild.audit <- '../../skyIslands_saved/data/checks/guild_id_review'
dir.create(dir.guild.audit, recursive=TRUE, showWarnings=FALSE)
source("src/guild_review.R")

## ***********************************************************************
## Parasite summaries at species, genus and site-round grains
## ***********************************************************************

parasite.site.sr.sp.sum <- spec.net %>%
    group_by(Site, Year, SampleRound, GenusSpecies, Genus) %>%
    summarise(SpAbundance = length(GenusSpecies),
              SpParasitism = sum.observed(ParasitePresence),
              SpCrithidiaPresence= sum.observed(CrithidiaPresence),
              SpApicystisSpp = sum.observed(ApicystisSpp),
              SpNosemaBombi= sum.observed(NosemaBombi),
              SpNosemaCeranae = sum.observed(NosemaCeranae),
              SpAscosphaeraSpp= sum.observed(AscosphaeraSpp),
              SpCrithidiaExpoeki = sum.observed(CrithidiaExpoeki),
              SpCrithidiaBombi = sum.observed(CrithidiaBombi),
              SpCrithidiaSpp = sum.observed(CrithidiaSpp),
              SpScreened = sum(Apidae == 1, na.rm=TRUE)
              )

parasite.site.sr.gen.sum <- spec.net %>%
    group_by(Site, Year, SampleRound, Genus) %>%
    summarise(GenusAbundance = length(GenusSpecies),
              GenusParasitism = sum.observed(ParasitePresence),
              GenusCrithidiaPresence= sum.observed(CrithidiaPresence),
              GenusApicystisSpp = sum.observed(ApicystisSpp),
              GenusNosemaBombi= sum.observed(NosemaBombi),
              GenusNosemaCeranae = sum.observed(NosemaCeranae),
              GenusAscosphaeraSpp= sum.observed(AscosphaeraSpp),
              GenusCrithidiaExpoeki = sum.observed(CrithidiaExpoeki),
              GenusCrithidiaBombi = sum.observed(CrithidiaBombi),
              GenusCrithidiaSpp = sum.observed(CrithidiaSpp),
              GenusScreened = sum(Apidae == 1, na.rm=TRUE)
              )

parasite.site.sr.sum <- spec.net %>%
    group_by(Site, Year, SampleRound) %>%
    summarise(SiteBeeAbundance = length(GenusSpecies[Family %in%
                                                     bee.families]),

              SiteBeeRichness= observed.richness(GenusSpecies[Family %in%
                                                          bee.families]),

              SiteBeeDiversity=observed.shannon(GenusSpecies[Family %in% bee.families]),

              SiteParasitism = sum.observed(ParasitePresence),
              SiteCrithidiaPresence= sum.observed(CrithidiaPresence),
              SiteApicystisSpp = sum.observed(ApicystisSpp),
              SiteNosemaBombi= sum.observed(NosemaBombi),
              SiteNosemaCeranae = sum.observed(NosemaCeranae),
              SiteAscosphaeraSpp= sum.observed(AscosphaeraSpp),
              SiteCrithidiaExpoeki = sum.observed(CrithidiaExpoeki),
              SiteCrithidiaBombi = sum.observed(CrithidiaBombi),
              SiteCrithidiaSpp = sum.observed(CrithidiaSpp),

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
## Site x Year x SampleRound plant-pollinator networks
## *******************************************************************

## Save plant x insect matrices for every guild, including verified empty
## rounds. Full network metrics are calculated in 8_network_metrics.R.
net.site.sr.guild.networks <- make_guild_networks(
  spec.net, net.site.sr.events %>% select(all_of(site.sr.keys)),
  guilds=pollinator_guilds(bee.families, butterfly.hawk.families),
  network_dir='../data/networks',
  metrics_dir='../data/splevel_network_metrics',
  audit_dir='../../skyIslands_saved/data/checks/network_id_review',
  compute_bee_metrics=FALSE)

## Retain only Chao1-estimated asymptotic plant partner richness per bee
## species and sampling round. Keep the legacy column name rare.degree and
## the previous calcSpec restriction to >= 3 taxa in both network guilds.
site.sr.sp.keys <- c(site.sr.keys, 'GenusSpecies')
net.site.sr.bee.partner.richness <- local({
  load('../data/networks/YearSR_PlantPollinator_Bees.Rdata')
  eligible <- which(vapply(nets, function(x) all(dim(x) > 2), logical(1)))

  ## One data frame per network; remove list names before binding rows.
  partner.richness <- lapply(unname(eligible), function(i) {
    data.frame(Site=unname(sites[i]), Year=as.numeric(years[i]),
               SampleRound=unname(SRs[i]),
               GenusSpecies=colnames(nets[[i]]),
               rare.degree=as.numeric(apply(nets[[i]], 2, fossil::chao1)))
  })

  ## Retain the join columns even if no network is eligible.
  empty <- data.frame(Site=character(), Year=numeric(), SampleRound=character(),
                      GenusSpecies=character(), rare.degree=numeric())
  bind_rows(c(list(empty), partner.richness))
})

spec.net <- left_join(spec.net, net.site.sr.bee.partner.richness,
                     by=site.sr.sp.keys, relationship='many-to-one')

## *******************************************************************
## Vegetation preparation
## *******************************************************************

veg <-
    read.csv("../../skyIslands_saved/data/relational/relational/traditional/veg-complete.csv",
             stringsAsFactors=FALSE, colClasses=c(SampleRound='character'))

veg$Date <- as.Date(veg$Date, format='%Y-%m-%d')
stopifnot(!anyNA(veg$Date))
veg$Year <- as.numeric(format(veg$Date, '%Y'))
veg$Quadrat <- trimws(as.character(veg$Quadrat))


## capture survey effort BEFORE filtering empty plant entries.
if (!'PlantGenusSpecies' %in% names(veg)) {
  veg$PlantGenusSpecies <- fix.white.space(paste(
    ifelse(is.na(veg$PlantGenus),'',veg$PlantGenus),
    ifelse(is.na(veg$PlantSpecies),'',veg$PlantSpecies),
    ifelse(is.na(veg$PlantVar),'',veg$PlantVar)))
}

## Vegetation was surveyed exactly once per Site x Year x SampleRound,
## including rounds without insect captures. Empty quadrats may be absent
## from the raw vegetation table. Reconstruct the fixed protocol below
## (20 quadrats in 2012; 27 later), not one set per netting condition.
veg.site.sr.quads.recorded <- unique(veg[, c('Site', 'Year', 'SampleRound', 'Quadrat'),
                             drop=FALSE])
## Retain unidentified flowering plants for abundance, but not richness.
veg$PlantGenusSpecies[is.na(veg$PlantGenusSpecies) |
                     !nzchar(trimws(veg$PlantGenusSpecies))] <- NA_character_


## All vegetation plant names were reconciled against the same
## GBIF mapping as insect visit plants in relational/1prep.R. 
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
if (any(invalid.quad)) warning('',
 sum(invalid.quad), ' floral rows have nonstandard quadrat IDs; ',
 'see checks/nonstandard_vegetation_quadrats.csv for review.')
veg <- veg[!invalid.quad, ,drop=FALSE]
## Each documented site-round had one vegetation survey. A round with
## no floral records is therefore a zero-flower survey, not missing effort.
veg.site.sr.rounds.recorded <- unique(veg.site.sr.quads.recorded[
  , c('Site', 'Year', 'SampleRound')])
## Do not turn vegetation without
## verified event geography into a newly inferred sampling location.
veg.rounds.without.net.events <- dplyr::anti_join(veg.site.sr.rounds.recorded,
 net.site.sr.events, by=site.sr.keys)
write.csv(veg.rounds.without.net.events,
 file='../../skyIslands_saved/data/checks/vegetation_without_verified_event.csv',
 row.names=FALSE)
if (nrow(veg.rounds.without.net.events)) warning('',
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
veg.site.sr.quads.recorded$RecordedVegetation <- TRUE
veg.site.sr.quads.design <- dplyr::left_join(veg.site.sr.quads.design,
  veg.site.sr.quads.recorded[, c('Site','Year','SampleRound','Quadrat',
                    'RecordedVegetation')],
  by=c('Site','Year','SampleRound','Quadrat'),
  relationship='one-to-one')
veg.site.sr.quads.design$RecordedVegetation[is.na(veg.site.sr.quads.design$RecordedVegetation)] <- FALSE
veg.site.sr.quads.design$ImputedZeroQuadrat <- !veg.site.sr.quads.design$RecordedVegetation

## Missing bloom counts on identified plant rows are not evidence of zero.
veg.quad.missing <- veg %>%
  filter(!is.na(PlantGenusSpecies), is.na(NumBlooms) | NumBlooms == '') %>%
  distinct(Site, Year, SampleRound, Quadrat) %>%
  mutate(MissingBloomCount=TRUE)

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
  summarise(FloralAbundance=sum(NumBloomsNum),
            FloweringPlantAbundance=sum(PlantCount), .groups='drop')

veg.site.sr.quad.sum <- veg.site.sr.quad.sp.sum %>%
  group_by(Quadrat, Site, Year, SampleRound) %>%
  summarise(FloralDiversity=shannon.counts(FloralAbundance[!is.na(PlantGenusSpecies)]),
            FloweringPlantDiversity=shannon.counts(FloweringPlantAbundance[!is.na(PlantGenusSpecies)]),
            FloralAbundance=sum(FloralAbundance),
            FloralRichness=observed.richness(PlantGenusSpecies),
            FloweringPlantAbundance=sum(FloweringPlantAbundance),
            .groups='drop')

## Compute Shannon before collapsing the abundance vectors (dplyr
## summarise can otherwise refer to the newly computed scalar).
veg.site.sr.quads.design <- left_join(
  veg.site.sr.quads.design, veg.site.sr.quad.sum %>% mutate(HasBloomRecords=TRUE),
  by=c(site.sr.keys, 'Quadrat'), relationship='one-to-one')

## All rounds were surveyed; absent bloom rows are zeros.
## Do not replace missing measurements within observed blooming quadrats.
empty.quad <- is.na(veg.site.sr.quads.design$HasBloomRecords)
veg.site.sr.quads.design[empty.quad, c('FloralAbundance', 'FloralRichness',
                                     'FloweringPlantAbundance')] <- 0
veg.site.sr.quads.design$HasBloomRecords <- NULL
veg.site.sr.quads.design <- left_join(veg.site.sr.quads.design, veg.quad.missing,
  by=c(site.sr.keys, 'Quadrat'), relationship='one-to-one')
missing.blooms <- veg.site.sr.quads.design$MissingBloomCount %in% TRUE
veg.site.sr.quads.design[missing.blooms, c('FloralAbundance', 'FloralRichness',
  'FloweringPlantAbundance', 'FloralDiversity', 'FloweringPlantDiversity')] <- NA
veg.site.sr.quads.design$MissingBloomCount <- NULL
write.csv(veg.site.sr.quads.design, file='../data/veg_quadrat_effort.csv', row.names=FALSE)

## Abundance/richness include every protocol quadrat. Shannon means are
## conditional on quadrats with calculable diversity.
veg.site.sr.sum <- veg.site.sr.quads.design %>%
    group_by(Site, Year, SampleRound) %>%
    summarise(MeanFloralAbundance = mean(FloralAbundance),
              MeanFloralRichness= mean(FloralRichness),
              MeanFloralDiversity =mean.observed(FloralDiversity),
              MeanFloweringPlantDiversity =mean.observed(FloweringPlantDiversity),
              MeanFloweringPlantAbundance =mean(FloweringPlantAbundance)
              )


write.csv(veg.site.sr.sum, file="../data/veg.csv", row.names=FALSE)

## Combine verified Net events, captures and independent vegetation summaries.
net.site.sr.sum <- net.site.sr.events %>%
  left_join(net.site.sr.sum, by=site.sr.keys, relationship='one-to-one') %>%
  left_join(veg.site.sr.sum, by=site.sr.keys, relationship='one-to-one')

## Only verified zero-capture rounds receive zeros for abundance/richness.
## Shannon diversity stays NA when there are no individuals.
## Parasite rates remain NA, and vegetation and conditions are unchanged.
empty.capture <- is.na(net.site.sr.sum$PollAbundance)
capture.cols <- grep(
  '^(Poll|Bee|Syrphid|HB|Bombus|NonBombusHB|ButterflyHawk|VisitedFloral).*(Abundance|Richness)$',
  names(net.site.sr.sum), value=TRUE)
for (col in capture.cols) net.site.sr.sum[[col]][empty.capture] <- 0
message('Documented Net rounds without specimens: ', sum(empty.capture))

write.csv(net.site.sr.sum, file='../data/sitestats.csv', row.names=FALSE)

## floral richness across the entire meadow across sampling rounds
veg.year.sum <- veg.blooming %>%
    group_by(Site, Year) %>%
    summarise(TotalFloralRichness=observed.richness(PlantGenusSpecies),
              HasBloomRecords=TRUE, .groups='drop')
veg.year.sum <- net.site.sr.events %>%
  distinct(Site, Year) %>%
  full_join(veg.year.sum, by=c('Site', 'Year'), relationship='one-to-one') %>%
  mutate(TotalFloralRichness=if_else(is.na(HasBloomRecords), 0L, TotalFloralRichness)) %>%
  select(-HasBloomRecords)
## A missing bloom measurement can prevent a complete annual richness estimate.
veg.year.sum <- veg.year.sum %>%
  left_join(veg.quad.missing %>% distinct(Site, Year) %>% mutate(MissingBloomCount=TRUE),
            by=c('Site', 'Year'), relationship='one-to-one') %>%
  mutate(TotalFloralRichness=if_else(MissingBloomCount %in% TRUE, NA_integer_, TotalFloralRichness)) %>%
  select(-MissingBloomCount)

write.csv(veg.year.sum, file="../data/veg_species_richness.csv",
          row.names=FALSE)

## *******************************************************************
## Merge site-round and species summaries onto Net specimens
## *******************************************************************

## Attach Site x Year x SampleRound summaries to Net specimen rows.
## Use a left join: no synthetic specimen rows should be created.
## Keep original specimen conditions under their original names; overlapping
## site-round summaries receive an explicit .site.sr suffix.
spec.net <- dplyr::left_join(spec.net, net.site.sr.sum, by=site.sr.keys,
                             suffix=c('', '.site.sr'),
                             relationship='many-to-one')

## Reuse species summaries calculated earlier, without a CSV round trip.
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
  message('', sum(missing.species),
          ' unidentified species-summary groups retained, ',
          'not matched as taxonomic species.')
}
net.site.sr.sp.stats.join <- net.site.sr.sp.stats[
  !missing.species, , drop=FALSE]

## Join identified species summaries without adding specimen rows.
spec.net <- dplyr::left_join(spec.net, net.site.sr.sp.stats.join, by=site.sr.sp.keys,
                             relationship='many-to-one')

spec.net$GenusSpecies[is.na(spec.net$GenusSpecies)] <- ""
save(spec.net, file="../data/spec_net.Rdata")

