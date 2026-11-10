setwd("../skyIslands_saved/data/relational/relational")
rm(list=ls())
library(RSQLite)

## connect to the relational database
con <- dbConnect(dbDriver("SQLite"), dbname='si.db')

## **************************************************
## specimen
## **************************************************

## make a table containing everything
##  Only verified conditions and geography are
## required for a specimen to enter specimens-complete.csv. An unknown
## plant or insect identification is NOT a reason to discard a specimen.
## Geography uses Site-level FKs, never subsite matching.
sql <- paste('SELECT * FROM tblSpecimens',
             'LEFT JOIN tblInsect',
             'ON tblSpecimens.InsectFK = tblInsect.InsectPK',
             'LEFT JOIN tblPlant',
             'ON tblSpecimens.PlantFK = tblPlant.PlantPK',
             'LEFT JOIN tblPan',
             'ON tblSpecimens.PanFK = tblPan.PanPK',
             'INNER JOIN tblConditions',
             'ON tblSpecimens.ConditionsFK = tblConditions.ConditionsPK',
             'INNER JOIN tblGeography',
             'ON tblSpecimens.GeographyFK = tblGeography.GeographyPK',
             'AND tblConditions.GeographyFK = tblGeography.GeographyPK')
res.complete <- dbGetQuery(con, sql)
names(res.complete) <- make.unique(names(res.complete))

##  Audit exclusions. Do not assign specimen localities
## when the sample event has no reliable conditions or site geography.
excluded.specimens <- dbGetQuery(con, paste(
  'SELECT s.* FROM tblSpecimens AS s',
  'LEFT JOIN tblConditions AS c ON s.ConditionsFK = c.ConditionsPK',
  'LEFT JOIN tblGeography AS g ON s.GeographyFK = g.GeographyPK',
  'AND c.GeographyFK = g.GeographyPK',
  'WHERE c.ConditionsPK IS NULL OR g.GeographyPK IS NULL'))
dir.create('../../checks/relational_audit', recursive=TRUE,
           showWarnings=FALSE)
write.csv(excluded.specimens,
  '../../checks/relational_audit/specimens_excluded_no_locality.csv',
  row.names=FALSE, na='')
original.n <- dbGetQuery(con,
  'SELECT COUNT(*) AS n FROM tblSpecimens')$n
if (nrow(res.complete) + nrow(excluded.specimens) != original.n)
  stop('Specimen exports do not reconcile with tblSpecimens. Check joins.')
message('Specimens excluded without trusted event/geography: ',
        nrow(excluded.specimens))

##  Event-anchored export preserves verified sampling
## events with no specimens. Such rows have geography and conditions,
## but NA specimen fields. This is NOT the specimen-counting table.
event.sql <- paste(
  'SELECT c.*, g.Site, g.Country, g.State, g.County, g.Meadow,',
  'g.Forest, g.MtRange, g.Lat, g.Long, g.Elev, g.Area,',
  's.*, i.*, p.* FROM tblConditions AS c',
  'INNER JOIN tblGeography AS g ON c.GeographyFK = g.GeographyPK',
  'LEFT JOIN tblSpecimens AS s ON s.ConditionsFK = c.ConditionsPK',
  'AND s.GeographyFK = c.GeographyFK',
  'LEFT JOIN tblInsect AS i ON s.InsectFK = i.InsectPK',
  'LEFT JOIN tblPlant AS p ON s.PlantFK = p.PlantPK')
event.specimens <- dbGetQuery(con, event.sql)
names(event.specimens) <- make.unique(names(event.specimens))
event.specimens$Year <- format(as.Date(event.specimens$Date, '%Y-%m-%d'),
                               '%Y')
event.specimens$HasSpecimen <- !is.na(event.specimens$UniqueID) &
  nzchar(trimws(as.character(event.specimens$UniqueID)))
write.csv(event.specimens,
  'traditional/conditions-specimens-complete.csv',
  row.names=FALSE, na='')

## drop unwanted columns
drop <- c('InsectPK', 'InsectFK',
          'GeographyPK', 'GeographyFK',
          'PlantPK', 'PlantFK',
          'PanPK', 'PanFK',
          'ConditionsPK', 'ConditionsFK')

##  Drop only columns that actually exist and
## remove provenance join fields without changing original site-level output.
res.complete <- res.complete[, !names(res.complete) %in% c(drop,
  'SourceDate', 'SourceSite', 'SourceSampleRound',
  'SourceMethod', 'SourceNetNumber'), drop=FALSE]

print(paste("after specimen traditional, dim=", nrow(res.complete)))

## Keep numeric missing values as NA in memory.
## Blank cells are written only in the CSV representation.

## write full table to csv
write.csv(res.complete, na='', file='traditional/specimens-complete.csv',
          row.names=FALSE)

## ## ## generate barcodes
## library(devtools)
## devtools::install_github("yihanwu/baRcodeR", build_vignettes = TRUE)
## library(baRcodeR)


## have.labels <- grepl("BBSL", res.complete$UniqueID) | as.Date(res.complete$Date, "%m/%d/%y") < "2013-01-01"

## custom_create_PDF(Labels = res.complete$UniqueID[!have.labels],
##                   name = "SI-labelsCustom", numcol=8)

## **************************************************
## veg
## **************************************************

## make a table containing everything
##  Unknown plant or survey keys remain visible
## in vegetation exports instead of silently dropping observations.
sqlveg <- paste('SELECT * FROM tblVeg',
             'LEFT JOIN tblQuad',
             'ON tblVeg.QuadFK = tblQuad.QuadPK',
             'LEFT JOIN tblPlant',
             'ON tblVeg.PlantFK = tblPlant.PlantPK',
             'LEFT JOIN tblSample',
             'ON tblVeg.SampleFK = tblSample.SamplePK'
             ## 'JOIN tblGeography',
             ## 'ON tblVeg.GeographyFK = tblGeography.GeographyPK'
             )
veg.complete <- dbGetQuery(con, sqlveg)
## UPDATED [CONS-20]: Lookup joins must preserve observation counts.
if (nrow(veg.complete) != dbGetQuery(con,
  'SELECT COUNT(*) AS n FROM tblVeg')$n)
  stop('Vegetation reconstruction changed the number of observations.')

## drop unwanted columns
drop <- c('QuadPK', 'QuadFK',
          'PlantPK', 'PlantFK',
          'SamplePK', 'SampleFK')

veg.complete <- veg.complete[, !names(veg.complete) %in% drop, drop=FALSE]

veg.complete <- veg.complete[,!grepl("\\..", colnames(veg.complete))]

print(paste("after veg traditional, dim=", nrow(veg.complete)))

## Preserve missing values as NA before CSV export.

## write full table to csv
write.csv(veg.complete, na='', file='traditional/veg-complete.csv',
          row.names=FALSE)



## **************************************************
## bloom
## **************************************************

## make a table containing everything
## Preserve unknown bloom plant/sample matches.
sqlbloom <- paste('SELECT * FROM tblBloom',
             'LEFT JOIN tblPlant',
             'ON tblBloom.PlantFK = tblPlant.PlantPK',
             'LEFT JOIN tblSample',
             'ON tblBloom.SampleFK = tblSample.SamplePK'
             ##   'JOIN tblGeography',
             ## 'ON tblBloom.GeographyFK = tblGeography.GeographyPK'
             )
bloom.complete <- dbGetQuery(con, sqlbloom)
## Lookup joins must preserve observation counts.
if (nrow(bloom.complete) != dbGetQuery(con,
  'SELECT COUNT(*) AS n FROM tblBloom')$n)
  stop('Bloom reconstruction changed the number of observations.')

## drop unwanted columns
drop <- c(
    'PlantPK', 'PlantFK',
    'SamplePK', 'SampleFK')

bloom.complete <- bloom.complete[, !names(bloom.complete) %in% drop, drop=FALSE]

bloom.complete <- bloom.complete[,!grepl("\\..", colnames(bloom.complete))]

print(paste("after bloom traditional, dim=", nrow(bloom.complete)))

##  Preserve missing values as NA before CSV export.

## write full table to csv
write.csv(bloom.complete, na='', file='traditional/bloom-complete.csv',
          row.names=FALSE)


## **************************************************
## close connection to database
dbDisconnect(con)
