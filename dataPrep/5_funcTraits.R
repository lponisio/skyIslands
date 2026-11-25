rm(list=ls())
library(FD)
library(readxl)
library(tidyverse)

setwd('~/')
source('lab_paths.R')
setwd(file.path(local.path, 'skyIslands/dataPrep'))
source('src/calcFuncUniqOrig.R')

load('../data/spec_net.Rdata')
site.sr.keys <- c('Site', 'Year', 'SampleRound')
sitestats <- read.csv('../data/sitestats.csv', colClasses=c(SampleRound='character'))
events <- sitestats %>% select(all_of(site.sr.keys)) %>% distinct()

bee.families <- c('Andrenidae', 'Apidae', 'Colletidae', 'Halictidae',
                  'Megachilidae', 'Melittidae')
spec.net <- spec.net %>% filter(Family %in% bee.families) %>%
  mutate(GenusSpecies=na_if(trimws(GenusSpecies), ''),
         Year=as.numeric(as.character(Year)), SampleRound=as.character(SampleRound))

## *******************************************************************
## Prepare species traits and one common trait space
## *******************************************************************

trait.path <- '../../PonisioLab/bee_traits/western_us_bee_traits.xlsx'
traits <- read_excel(trait.path, sheet='Data') %>%
  mutate(across(where(is.character), ~na_if(trimws(.x), '')))

##  Each trait block (nesting, size, pollen carrying, reproductive
## strategy/sociality) has weight 1.

bee.traits <- c('NestLocation', 'NestPartitions', 'NestConstruction',
                'NestLinning', 'PrimaryNestMaterial', 'MeanITD',
                'PollenCarry', 'ReproStrat', 'Sociality')
bee.weights <- setNames(c(rep(1/5, 5), 1, 1, 1/2, 1/2), bee.traits)

## Review duplicate records and trait coverage.
audit.dir <- '../../skyIslands_saved/data/checks/functional_traits'
source('src/trait_review.R')

## Species x trait matrix for the observed bee species.
## Gower distance handles partially missing traits pairwise.
trait.space <- traits %>%
  filter(GenusSpecies %in% trait.audit$GenusSpecies[trait.audit$InTraitSpace]) %>%
  select(GenusSpecies, all_of(bee.traits)) %>%
  mutate(across(-c(GenusSpecies, MeanITD), as.factor)) %>%
  column_to_rownames('GenusSpecies')

weights <- bee.weights

## *******************************************************************
## Community matrices: Site, Site x Year, Site x Year x SampleRound
## *******************************************************************

scales <- list(Site='Site', SiteYear=c('Site', 'Year'),
               SiteYearSR=site.sr.keys)
comms <- lapply(names(scales), function(suffix) {
  make.comm.matrix(spec.net, events, scales[[suffix]],
                   rownames(trait.space), suffix,
                   coverage=(suffix == 'SiteYearSR'))
})
names(comms) <- names(scales)
a <- do.call(rbind, lapply(comms, function(x) x$a))

## One regional trait space for species roles at all three scales.
## Calculate community FD only for site-year-sampling rounds.
out <- calcFuncUniqOrigComm(trait.space, weights, a, w.abun=TRUE,
                           corr='cailliez',
                           fd.ids=rownames(comms$SiteYearSR$a))

## *******************************************************************
## Export community and species-level metrics and attach them to specimens
## *******************************************************************

spec.net <- spec.net %>%
  unite('SiteYearSR', Site, Year, SampleRound, sep=';', remove=FALSE) %>%
  unite('SiteYear', Site, Year, sep=';', remove=FALSE) %>%
  left_join(traits, by='GenusSpecies', suffix=c('', '.trait'),
            na_matches='never', relationship='many-to-one')

## Attach species functional roles at all three scales.
for (suffix in names(scales)) {
  keys <- scales[[suffix]]
  species <- comms[[suffix]]$groups %>%
    select(all_of(keys), commID) %>%
    inner_join(out$by.comm.mets, by='commID', relationship='one-to-many') %>%
    select(-commID) %>%
    rename_with(~paste0(.x, suffix), c(uniq, originality, unweightedOriginality))
  spec.net <- spec.net %>%
    left_join(species, by=c(keys, 'GenusSpecies'), relationship='many-to-one')
}

## Append community FD and trait coverage to sitestats, retaining all rounds.
fd.cols <- c('BeeFDis', 'BeeFEve', 'BeeFDiv',
             'BeeTraitAbundance', 'BeeTraitRichness', 'BeeTraitCoverage')
site.round.fd <- comms$SiteYearSR$groups %>%
  left_join(out$fd, by='commID', relationship='one-to-one') %>%
  select(all_of(c(site.sr.keys, fd.cols)))

## Replace existing FD columns if the script is rerun.
sitestats <- sitestats %>%
  select(-any_of(fd.cols)) %>%
  left_join(site.round.fd, by=site.sr.keys, relationship='one-to-one')
write.csv(sitestats, '../data/sitestats.csv', row.names=FALSE)

## Retain site-round FD columns on specimens for existing analyses.
spec.net <- spec.net %>%
  left_join(site.round.fd %>%
              rename_with(~paste0(.x, 'SiteYearSR'), all_of(fd.cols)),
            by=site.sr.keys, relationship='many-to-one')

spec.net$GenusSpecies[is.na(spec.net$GenusSpecies)] <- ''
save(spec.net, file='../data/spec_traits.Rdata')
write.csv(spec.net, '../data/spec_traits.csv', row.names=FALSE)
