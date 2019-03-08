setwd("~/Dropbox/skyIslands_saved/data/relational/relational")
rm(list=ls())
library(RSQLite)

## connect to the relational database
con <- dbConnect(dbDriver("SQLite"), dbname='si.db')

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

print(paste("after traditional, dim=", nrow(res.complete)))

## set NA to blanks
res.complete[is.na(res.complete)] <- ''

## write full table to csv
write.csv(res.complete, file='traditional/specimens-complete.csv',
          row.names=FALSE)

## ## ## generate barcodes
library(devtools)
## devtools::install_github("yihanwu/baRcodeR", build_vignettes = TRUE)
## library(baRcodeR)


## have.labels <- grepl("BBSL", res.complete$UniqueID) | as.Date(res.complete$Date, "%m/%d/%y") < "2013-01-01"

## custom_create_PDF(Labels = res.complete$UniqueID[!have.labels],
##                   name = "SI-labelsCustom", numcol=8)

## **************************************************
## close connection to database
dbDisconnect(con)
