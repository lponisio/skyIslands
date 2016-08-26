setwd("~/Dropbox/SkyIslands/data/relational/data/relational")
rm(list=ls())
library(RSQLite)

## connect to the relational database
con <- dbConnect(dbDriver("SQLite"), dbname='si.db')

## **************************************************
## make a table containing everything
sql <- paste("SELECT * FROM tblSpecimens",
             "JOIN tblInsect",
             "ON tblSpecimens.InsectFK = tblInsect.InsectPK",
             "JOIN tblPlant",
             "ON tblSpecimens.PlantFK = tblPlant.PlantPK",
             "JOIN tblConditions",
             "ON tblSpecimens.ConditionsFK = tblConditions.ConditionsPK",
             "JOIN tblGeography",
             "ON tblSpecimens.GeographyFK = tblGeography.GeographyPK")
res.complete <- dbGetQuery(con, sql)

## check dimensions of dataset
dim(res.complete)

## drop unwanted columns
drop <- c("InsectPK", "InsectFK",
          "GeographyPK", "GeographyFK",
          "PlantPK", "PlantFK",
          "ConditionsPK", "ConditionsFK")
res.complete <- res.complete[-match(drop, names(res.complete))]
drop <- c("GeographyFK")
res.complete <- res.complete[-match(drop, names(res.complete))]

## set NA to blanks
res.complete[is.na(res.complete)] <- ""

## write full table to csv
write.csv(res.complete, file="traditional/specimens-complete.csv",
          row.names=FALSE)
## **************************************************
## close connection to database
sqliteCloseConnection(con)
