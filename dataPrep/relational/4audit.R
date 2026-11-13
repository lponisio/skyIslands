
## *******************************************************
## optional relational integrity audit (run after 3join.R)
## *******************************************************
library(RSQLite)
setwd('../skyIslands_saved/data/relational/relational')
audit.dir <- '../../checks/relational_audit'
dir.create(audit.dir, recursive=TRUE, showWarnings=FALSE)
con <- dbConnect(dbDriver('SQLite'), dbname='si.db')

## Count records at each stage, including the traditional exports.
record_count <- function(table) as.integer(dbGetQuery(con,
  paste0('SELECT COUNT(*) AS n FROM ', table))$n)
tables <- c('tblSpecimens', 'tblConditions', 'tblGeography',
            'tblInsect', 'tblPlant', 'tblPan', 'tblVeg',
            'tblBloom', 'tblSample', 'tblQuad')
counts <- data.frame(Table=tables,
                     N=vapply(tables, record_count, integer(1)))
write.csv(counts, file=file.path(audit.dir, 'table_row_counts.csv'),
          row.names=FALSE)

##  Missing specimen conditions/geography are a
## deliberate exclusion from locality-labelled specimen exports.
## Unmatched plant IDs remain descriptive and do not trigger exclusion.
spec <- dbReadTable(con, 'tblSpecimens')
fks <- c('InsectFK', 'PlantFK', 'PanFK',
         'ConditionsFK', 'GeographyFK')
fk.audit <- data.frame(Field=fks,
                       Missing=vapply(fks, function(x) sum(is.na(spec[[x]])),
                                      integer(1)),
                       Total=nrow(spec))
write.csv(fk.audit, file=file.path(audit.dir, 'missing_specimen_keys.csv'),
          row.names=FALSE)
for (fk in fks) {
  bad <- spec[is.na(spec[[fk]]),
              intersect(c('UniqueID', 'TempID', 'SourceSite',
                          'SourceDate', 'SourceMethod', 'SourcePlantName',
                          fk), names(spec)), drop=FALSE]
  write.csv(bad, file=file.path(audit.dir, paste0('missing_', fk, '.csv')),
            row.names=FALSE)
}

## Report duplicate IDs, rather than silently removing one of them.
dups <- duplicated(spec$UniqueID) |
        duplicated(spec$UniqueID, fromLast=TRUE)
write.csv(spec[dups, ], file=file.path(audit.dir,
          'duplicate_specimen_IDs.csv'), row.names=FALSE)

## Additional foreign keys for plant/quadrat observation fact tables.
for (table in c('tblVeg', 'tblBloom')) {
  dat <- dbReadTable(con, table)
  for (fk in intersect(c('QuadFK', 'PlantFK', 'SampleFK'), names(dat))) {
    out <- dat[is.na(dat[[fk]]), , drop=FALSE]
    write.csv(out, file=file.path(audit.dir,
      paste0(table, '_missing_', fk, '.csv')), row.names=FALSE)
  }
}

##  The specimen export intentionally excludes records
## without matched conditions/geography. Reconcile source, kept and dropped.
exports <- c('specimens-complete.csv', 'veg-complete.csv',
             'bloom-complete.csv')
observed <- vapply(exports, function(f) {
  x <- read.csv(file.path('traditional', f), stringsAsFactors=FALSE)
  nrow(x)
}, integer(1))
excluded <- read.csv(file.path(audit.dir,
  'specimens_excluded_no_locality.csv'), stringsAsFactors=FALSE)
truth <- c(record_count('tblSpecimens') - nrow(excluded),
           record_count('tblVeg'), record_count('tblBloom'))
reconstruction <- data.frame(File=exports, Expected=truth,
                              Exported=observed, Preserved=(truth==observed))
write.csv(reconstruction, file=file.path(audit.dir,
          'reconstruction_counts.csv'), row.names=FALSE)
if (any(!reconstruction$Preserved))
  stop('Traditional exports unexpectedly changed row counts.')
##  Check every condition event is represented by one
## or more rows, including an NA-specimen row for an empty event.
ledger <- read.csv('traditional/conditions-specimens-complete.csv',
                   stringsAsFactors=FALSE)
expected.events <- dbGetQuery(con, paste(
  'SELECT COUNT(*) AS n FROM tblConditions AS c',
  'INNER JOIN tblGeography AS g ON c.GeographyFK = g.GeographyPK'))$n
ledger.events <- length(unique(ledger$ConditionsPK))
if (ledger.events != expected.events)
  stop('Condition-led export lost one or more geographic events.')
if (any(ledger$HasSpecimen == FALSE & !is.na(ledger$UniqueID)))
  stop('Empty sampling-event rows have a nonmissing UniqueID.')

dbDisconnect(con)
message('## Relational audit written to: ', normalizePath(audit.dir))
