## *******************************************************************
## Bee trait review
## Sourced by 5_funcTraits.R after loading traits and defining bee.traits.
## Uses traits, spec.net, bee.traits and audit.dir.
## *******************************************************************

dir.create(audit.dir, recursive=TRUE, showWarnings=FALSE)

## Prefer Sky Islands records, then the most complete trait record.
## Keep the original duplicate records for review.
traits <- traits %>%
  mutate(nTraits=rowSums(!is.na(pick(all_of(bee.traits)))),
         skyIslands=grepl('(^|, *)SI( *,|$)', coalesce(Project, ''))) %>%
  arrange(GenusSpecies, desc(skyIslands), desc(nTraits))

write.csv(traits %>% add_count(GenusSpecies) %>% filter(n > 1),
          file.path(audit.dir, 'duplicate_trait_records.csv'), row.names=FALSE)

traits <- traits %>% distinct(GenusSpecies, .keep_all=TRUE)

## Count available traits for each captured bee taxon.
trait.audit <- spec.net %>%
  count(GenusSpecies, name='BeeAbundance') %>%
  left_join(traits, by='GenusSpecies', na_matches='never',
            relationship='many-to-one') %>%
  mutate(BeeTraitsAvailable=rowSums(!is.na(pick(all_of(bee.traits)))),
         InTraitSpace=!is.na(GenusSpecies) & BeeTraitsAvailable > 0)

write.csv(trait.audit,
          file.path(audit.dir, 'trait_coverage_by_taxon.csv'), row.names=FALSE)

## Unmatched named taxa and unidentified individuals remain in spec.net,
## but cannot enter the functional trait space without trait information.
missing.taxa <- trait.audit %>%
  filter(!is.na(GenusSpecies), BeeTraitsAvailable == 0)

if (nrow(missing.taxa))
  warning('Named bee taxa without usable traits: ',
          paste(missing.taxa$GenusSpecies, collapse=', '),
          '. See trait_coverage_by_taxon.csv.')

## Check the observed species pool rather than unrelated workbook taxa.
## Do not silently remove trait columns or change their weights.
observed.traits <- traits %>%
  filter(GenusSpecies %in% trait.audit$GenusSpecies[trait.audit$InTraitSpace]) %>%
  select(all_of(bee.traits))

if (!nrow(observed.traits))
  stop('No captured bee taxa have usable traits.')

if (!is.numeric(observed.traits$MeanITD))
  stop('MeanITD must be numeric.')

n.values <- vapply(observed.traits, function(x)
  length(unique(na.omit(x))), integer(1))

if (any(n.values == 0) || n.values['MeanITD'] < 2)
  stop('A trait is entirely missing or MeanITD has no variation in the ',
       'observed trait pool. Review the trait workbook.')

rm(observed.traits, missing.taxa, n.values)
