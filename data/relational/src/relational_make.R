## *******************************************************
## create relational database
## *******************************************************
rm(list=ls())
setwd("~/Dropbox/SkyIslands/data/relational/data/relational")
library(RSQLite)
source('~/Dropbox/SkyIslands/analysis/data/src/misc.R')

conditions <- read.csv("../original/weather.csv", as.is=TRUE)
specimens <- read.csv("../original/specimens.csv", as.is=TRUE)
Geography <- read.csv("../original/geography.csv", as.is=TRUE)

## *******************************************************
## Start by importing the conditions information, from the conditions
## file.
## *******************************************************
## check that if there is already a database, remove it
if(file.exists("si.db")) file.remove("si.db") 
con <- dbConnect(dbDriver("SQLite"), dbname='si.db') 

## *******************************************************
## 1. Geographic infomation
## *******************************************************

keep <- c("Site", "Country", "State", "County", "Meadow", "Forest",
          "MtRange", "Lat", "Long")
geography <- unique(Geography[keep])
## next sort into alphabetical order
geography <- geography[match(sort(geography$Site), geography$Site),]

## generate primary geography key
geography <- cbind(GeographyPK=seq_len(nrow(geography)), geography)
rownames(geography) <- NULL
dbWriteTable(con, "tblGeography", geography, row.names=FALSE)

## Propagate geography key to the conditions table.
conditions$GeographyFK <-
  geography$GeographyPK[match(conditions$Site, geography$Site)]
## Propagate geography key to the specimens table.
specimens$GeographyFK <- geography$GeographyPK[match(specimens$Site,
                                                     geography$Site)] 

## write a .csv version of this table (just for ease of viewing)
write.csv(dbReadTable(con, "tblGeography"),
          file="tables/geography.csv", row.names=FALSE)

dbListTables(con)
## *******************************************************
## 2. Conditions
## *******************************************************

## Temporarily identify unique combinations:
keep <- c("GeographyFK", "Date")
conditions$cond.code <- apply(conditions[keep], 1, paste, collapse=";")
specimens$cond.code <- apply(specimens[keep], 1, paste, collapse=";")

## make table
keep <- c("Date", "Dos", "SampleRound", "StartTime", "EndTime",
          "TempStart", "TempEnd", "WindStart", "WindEnd", "SkyStart",
          "SkyEnd","GeographyFK", "cond.code")

cond <- unique(conditions[keep])
rownames(cond) <- NULL
cond <- cbind(ConditionsPK=seq_len(nrow(cond)), cond)
## Don't upload the cond.code column
dbWriteTable(con, "tblConditions", cond[-ncol(cond)], row.names=FALSE)

## Propagate conditions key to the conditions table.
conditions$ConditionsFK <-
  cond$ConditionsPK[match(conditions$cond.code, cond$cond.code)]
## Propagate conditions key to the specimens table.
specimens$ConditionsFK <-
  cond$ConditionsPK[match(specimens$cond.code, cond$cond.code)]

write.csv(dbReadTable(con, "tblConditions"),
          file="tables/conditions.csv", row.names=FALSE)

## *******************************************************
## 3. Insect species:
## *******************************************************

keep <- c("Order", "Family", "Genus", "SubGenus", "Species",
          "SubSpecies", "Determiner", "Author", "DateDetermined")

insects <- specimens[keep]
insects <- unique(insects)
insects$gen.sp <- paste(insects$Order,
                        insects$Family,
                        insects$Genus,
                        insects$SubGenus,
                        insects$Species,
                        insects$SubSpecies, sep=";") 
insects <- cbind(InsectPK=seq_len(nrow(insects)), insects)
rownames(insects) <- NULL

dbWriteTable(con, "tblInsect", insects[-ncol(insects)],
             row.names=FALSE)

## Propagate insect key to the specimens table.
specimens$gen.sp <- paste(specimens$Order,
                         specimens$Family,
                         specimens$Genus,
                         specimens$SubGenus,
                         specimens$Species,
                         specimens$SubSpecies, sep=";") 
specimens$InsectFK <- insects$InsectPK[match(specimens$gen.sp,
                                            insects$gen.sp)]

write.csv(dbReadTable(con, "tblInsect"),
          file="tables/insect.csv", row.names=FALSE)

## *******************************************************
## 4. Plant species:
## *******************************************************

keep <- c("FinalPlantSp", "PlantFamily") 
plants <- specimens[keep]
plants <- unique(plants)

PlantFamily <-plants$PlantFamily
PlantGenus <- sapply(strsplit(plants[,1], " "), function(x) x[1])
PlantGenus[is.na(PlantGenus)] <- ""
PlantSpecies <- sapply(strsplit(plants[,1], " "), function(x) x[2])
PlantSpecies[is.na(PlantSpecies)] <- ""
PlantVar <- sapply(strsplit(plants[,1], " "), function(x) x[3])
PlantVar[is.na(PlantVar)] <- ""
PlantSubSpecies <- sapply(strsplit(plants[,1], " "), function(x) x[4])
PlantSubSpecies[is.na(PlantSubSpecies)] <- ""


plants <- data.frame(PlantPK=seq_along(PlantGenus),
                     PlantFamily, PlantGenus, PlantSpecies,
                     PlantVar, PlantSubSpecies)

rownames(plants) <- NULL
dbWriteTable(con, "tblPlant", plants, row.names=FALSE)

## Propagate plant key to the specimens table.

specimens.plant.sp <- specimens[keep][[1]]
plants.plant.sp <- fix.white.space(paste(plants$PlantGenus,
                                         plants$PlantSpecies, plants$PlantVar,
                         plants$PlantSubSpecies))

specimens$PlantFK <- plants$PlantPK[match(specimens.plant.sp,
                                         plants.plant.sp)]

write.csv(dbReadTable(con, "tblPlant"),
          file="tables/plant.csv", row.names=FALSE)


## 6. Specimens:
## *******************************************************

keep <- c( "UniqueID",
          "TempID",
          "Collector",
          "SampleRound",
          "InsectFK",
          "PlantFK",
          "ConditionsFK",
          "GeographyFK")

specimens <- specimens[keep]
specimens <- unique(specimens)
rownames(specimens) <- NULL

dbWriteTable(con, "tblSpecimens", specimens, row.names=FALSE)

write.csv(dbReadTable(con, "tblSpecimens"),
          file="tables/specimens.csv", row.names=FALSE)

## *******************************************************
## close connection to database
## *******************************************************
sqliteCloseConnection(con)
