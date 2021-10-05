setwd("../skyIslands_saved/data/relational/relational")
rm(list=ls())
library(RSQLite)

## connect to the relational database
con <- dbConnect(dbDriver("SQLite"), dbname='si.db')

## **************************************************
## specimen
## **************************************************

## make a table containing everything
sql <- paste('SELECT * FROM tblSpecimens',
             'JOIN tblInsect',
             'ON tblSpecimens.InsectFK = tblInsect.InsectPK',
             'JOIN tblPlant',
             'ON tblSpecimens.PlantFK = tblPlant.PlantPK',
             'JOIN tblPan',
             'ON tblSpecimens.PanFK = tblPan.PanPK',
             'JOIN tblConditions',
             'ON tblSpecimens.ConditionsFK = tblConditions.ConditionsPK',
             'JOIN tblGeography',
             'ON tblSpecimens.GeographyFK = tblGeography.GeographyPK')
res.complete <- dbGetQuery(con, sql)

## drop unwanted columns
drop <- c('InsectPK', 'InsectFK',
          'GeographyPK', 'GeographyFK',
          'PlantPK', 'PlantFK',
          'PanPK', 'PanFK',
          'ConditionsPK', 'ConditionsFK')

res.complete <- res.complete[-match(drop, names(res.complete))]

res.complete <- res.complete[,!grepl("\\..", colnames(res.complete))]

print(paste("after specimen traditional, dim=", nrow(res.complete)))

## set NA to blanks
res.complete[is.na(res.complete)] <- ''

## write full table to csv
write.csv(res.complete, file='traditional/specimens-complete.csv',
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
sqlveg <- paste('SELECT * FROM tblVeg',
             'JOIN tblQuad',
             'ON tblVeg.QuadFK = tblQuad.QuadPK',
             'JOIN tblPlant',
             'ON tblVeg.PlantFK = tblPlant.PlantPK',
             'JOIN tblSample',
             'ON tblVeg.SampleFK = tblSample.SamplePK'
             ## 'JOIN tblGeography',
             ## 'ON tblVeg.GeographyFK = tblGeography.GeographyPK'
             )
veg.complete <- dbGetQuery(con, sqlveg)

## drop unwanted columns
drop <- c('QuadPK', 'QuadFK',
          'PlantPK', 'PlantFK',
          'SamplePK', 'SampleFK')

veg.complete <- veg.complete[-match(drop, names(veg.complete))]

veg.complete <- veg.complete[,!grepl("\\..", colnames(veg.complete))]

print(paste("after veg traditional, dim=", nrow(veg.complete)))

## set NA to blanks
veg.complete[is.na(veg.complete)] <- ''

## write full table to csv
write.csv(veg.complete, file='traditional/veg-complete.csv',
          row.names=FALSE)



## **************************************************
## bloom
## **************************************************

## make a table containing everything
sqlbloom <- paste('SELECT * FROM tblBloom',
             'JOIN tblPlant',
             'ON tblBloom.PlantFK = tblPlant.PlantPK',
             'JOIN tblSample',
             'ON tblBloom.SampleFK = tblSample.SamplePK'
             ##   'JOIN tblGeography',
             ## 'ON tblBloom.GeographyFK = tblGeography.GeographyPK'
             )
bloom.complete <- dbGetQuery(con, sqlbloom)

## drop unwanted columns
drop <- c(
    'PlantPK', 'PlantFK',
    'SamplePK', 'SampleFK')

bloom.complete <- bloom.complete[-match(drop, names(bloom.complete))]

bloom.complete <- bloom.complete[,!grepl("\\..", colnames(bloom.complete))]

print(paste("after bloom traditional, dim=", nrow(bloom.complete)))

## set NA to blanks
bloom.complete[is.na(bloom.complete)] <- ''

## write full table to csv
write.csv(bloom.complete, file='traditional/bloom-complete.csv',
          row.names=FALSE)


## **************************************************
## close connection to database
dbDisconnect(con)
