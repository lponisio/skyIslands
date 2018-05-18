## **********************************
## relational_prep.R script ##
## **********************************
## *******************************************************
## make files for use in relational database
## *******************************************************
rm(list=ls())
setwd("~/Dropbox/SkyIslands/data")

## load data
spec.data <- read.csv("raw/specimens.csv",  stringsAsFactors=FALSE,
                      colClasses=c("SampleRound"="character"))

## raw data was entered eith different names than in 2012 but they are
## what we want for the relational database

## *******************************************************
## next we will create the data structures that that we
## will use to construct our relational database
## *******************************************************

#spec.data <- spec.data[rep(seq_len(nrow(spec.data)), spec.data$SpecimenCount),]

#spec.data$SpecimenCount <- NULL
spec.data$SiteSubSite <-  spec.data$Site
spec.data$Site <- gsub("[0-9]", "", spec.data$Site)

print(paste("original number of specimens", nrow(spec.data)))

## add subsite one to all 2012 data
spec.2012 <- as.Date(spec.data$Date, "%m/%d/%y") < "2013-01-01"
spec.2017 <- as.Date(spec.data$Date, "%m/%d/%y") > "2013-01-01"

spec.data$SiteSubSite[spec.2012] <- paste0(spec.data$SiteSubSite[spec.2012], "1")

plant.keys <- read.csv("raw/plants.csv", stringsAsFactors=FALSE)

spec.data$FinalPlantSp[spec.2017] <- plant.keys$GenusSpecies[match(spec.data$FieldPlantID[spec.2017],
                                                                   plant.keys$FieldPlantID)]

spec.data$UniqueID[spec.2017] <- spec.data$TempID[spec.2017] <-
  1:length(spec.data$UniqueID[spec.2017])

spec.data$FinalPlantSp[is.na(spec.data$FinalPlantSp)] <- ""

spec.data$NetNumber <-
  sapply(strsplit(spec.data$SampleRound, ".", fixed=TRUE), function(x) x[2])

spec.data$SampleRound <-
  sapply(strsplit(spec.data$SampleRound, ".", fixed=TRUE), function(x) x[1])

spec.data$Method <- "Net"
spec.data$Method[is.na(spec.data$FieldPlantID)] <- "Pan"
spec.data$Method[spec.data$SampleRound == 0] <- "Vane"

check.data.spec <- aggregate(spec.data$Date, list(site = spec.data$Site,
                                                  date =
                                                    spec.data$Date),
                             length)


write.csv(spec.data, file="relational/data/original/specimens.csv",
          row.names=FALSE)

#source('../dataEntry/speciesIDs/AssignSpecies.R')

#write.csv(spec.data, file="relational/data/original/specimens.csv",
#          row.names=FALSE)


## *******************************************************
## create conditions file

data.weather <- read.csv("raw/weather.csv")

data.weather$SiteSubSite <-  paste0(data.weather$Site,
                                    data.weather$SubSite)

check.data.weather <- aggregate(data.weather$StartTime,
                                list(site = data.weather$Site,
                                     day = data.weather$SampleRound,
                                     date = data.weather$Date), length)

## write unique data to a table
write.csv(unique(data.weather), file="relational/data/original/weather.csv",
          row.names=FALSE)
## *******************************************************


geo <- read.csv("raw/geography.csv")

## data.geo <- data.frame(Site=geo$site,
##                    Meadow=geo$meadow,
##                    MtRange=geo$MtRange,
##                    Forest=geo$forest,
##                    County=geo$county,
##                    State=geo$state,
##                    Country=geo$country,
##                    Lat=geo$lat,
##                    Long=geo$long)

colnames(geo) <- gsub("\\.", "", colnames(geo))

geo$Site <- gsub("[0-9]", "", geo$Location_Name)

geo$SubSite <- substr(geo$Location_Name, 3, 3)

highergeog <- data.frame(Site = c("JC", "VC", "JM", "SC", "MM", "SM",
                                  "CC", "UK", "SS", "CH", "RP", "PL",
                                  "HM"),
                         State = c("NM", "NM", "NM", "NM", "NM", "NM",
                                   "NM", "NM", "NM", "AZ", "AZ", "AZ",
                                   "AZ"),
                         County = c("San Miguel", "Sandoval",
                                    "Rio Arriba", "Bernalillo",
                                    "Socorro", "Cibola", "Otero",
                                    "Otero", "Otero", "Cochise",
                                    "Cochise", "Graham", "Greenlee"),
                         Meadow = c("Jack's Creek Campground",
                                    "Valles Caldera National Preserve",
                                    "Cerro Pelon", "Kiwanis Meadow",
                                    "South Baldy Meadow",
                                    "La Mosca Lookout Tower Meadow",
                                    "Cathey Canyon",
                                    "Upper Karr Campground",
                                    "2 mi E of Sunspot",
                                    "Barfoot Park", "Rustler Park",
                                    "Hospital Flat Campground",
                                    "Hannagan Meadow"),
                         Forest = c("Santa Fe",
                                    "Valles Caldera National Preserve",
                                    "Santa Fe", "Cibola", "Cibola",
                                    "Cibola", "Lincoln", "Lincoln",
                                    "Lincoln", "Coronado", "Coronado",
                                    "Coronado", "Apache-Sitgreaves"),
                         MtRange=c("Pecos", "Jemez", "Jemez", "Sandias",
                                   "Magdalena", "San Mateo", "Sacramento",
                                   "Sacramento", "Sacramento",
                                   "Chiricahua", "Chiricahua",
                                   "Pina Leno", "Gila"))

geo <- merge(geo, highergeog, key="Site", all.x=TRUE)

data.geo <- data.frame(Site=geo$Site,
                       SubSite=geo$SubSite,
                       SiteSubSite=geo$Location_Name,
                       MtRange=geo$MtRange,
                       Forest=geo$Forest,
                       Meadow=geo$Meadow,
                       County=geo$County,
                       State=geo$State,
                       Country=geo$Country,
                       Lat=geo$DecimalLat,
                       Long=geo$DecimalLon,
                       Elev=geo$Elev0)

## write unique data to a table
write.csv(data.geo, file="relational/data/original/geography.csv",
          row.names=FALSE)

## **********************************
## relational_make.R script ##
## **********************************
  ## *******************************************************
  ## create relational database
  ## *******************************************************
  rm(list=ls())
setwd("~/Dropbox/SkyIslands/data/relational/data/relational")
library(RSQLite)
source('~/Dropbox/SkyIslands/analysis/data/src/misc.R')

conditions <- read.csv("../original/weather.csv", as.is=TRUE)
specimens <- read.csv("../original/specimens.csv", as.is=TRUE)
geo <- read.csv("../original/geography.csv", as.is=TRUE)

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

keep <- c("SiteSubSite", "Site", "Country", "State", "County",
          "Meadow", "Forest",
          "MtRange", "Lat", "Long")

geography <- unique(geo[keep])
## next sort into alphabetical order
geography <- geography[match(sort(geography$SiteSubSite), geography$SiteSubSite),]

## generate primary geography key
geography <- cbind(GeographyPK=seq_len(nrow(geography)), geography)
rownames(geography) <- NULL
dbWriteTable(con, "tblGeography", geography, row.names=FALSE)

## Propagate geography key to the conditions table.
conditions$GeographyFK <-
  geography$GeographyPK[match(conditions$SiteSubSite,
                              geography$SiteSubSite)]

## Propagate geography key to the specimens table.
specimens$GeographyFK <- geography$GeographyPK[match(specimens$SiteSubSite,
                                                     geography$SiteSubSite)]

## write a .csv version of this table (just for ease of viewing)
write.csv(dbReadTable(con, "tblGeography"),
          file="tables/geography.csv", row.names=FALSE)

dbListTables(con)

## *******************************************************
## 2. Conditions
## *******************************************************

## Temporarily identify unique combinations:
keep <- c("GeographyFK", "Date", "Method", "NetNumber", "Collector")
conditions$cond.code <- apply(conditions[keep], 1, paste, collapse=";")
specimens$cond.code <- apply(specimens[keep], 1, paste, collapse=";")

## make table
keep <- c("Date", "Dos", "Collector", "SampleRound", "NetNumber", "Method", "StartTime", "EndTime",
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

specimens$UniqueID[is.na(specimens$ConditionsFK)]

## specimens[specimens$UniqueID == "9296",]
## conditions[conditions$Date == "7/25/17",]

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

specimens$UniqueID[is.na(specimens$InsectFK)]

write.csv(dbReadTable(con, "tblInsect"),
          file="tables/insect.csv", row.names=FALSE)

## *******************************************************
## 4. Plant species:
## *******************************************************

keep <- c('FinalPlantSp')
plants <- specimens[keep][[1]]
plants <- sort(unique(plants))

PlantGenus <- sapply(strsplit(plants, ' '), function(x) x[1])
PlantGenus[is.na(PlantGenus)] <- ''

PlantSpecies <- sapply(strsplit(plants, ' '), function(x) x[2])
PlantSpecies[is.na(PlantSpecies)] <- ''

PlantVar <- sapply(strsplit(plants, ' '), function(x) x[3])
PlantVar[is.na(PlantVar)] <- ''

PlantSubSpecies <- sapply(strsplit(plants, ' '), function(x) x[4])
PlantSubSpecies[is.na(PlantSubSpecies)] <- ''

plants <- data.frame(PlantPK=seq_along(PlantGenus),
                     PlantGenus, PlantSpecies, PlantVar,
                     PlantSubSpecies)

rownames(plants) <- NULL
dbWriteTable(con, 'tblPlant', plants, row.names=FALSE)

## Propagate plant key to the specimens table.
specimens.plant.sp <- specimens[keep][[1]]

plants.plant.sp <- fix.white.space(paste(plants$PlantGenus,
                                         plants$PlantSpecies,
                                         plants$PlantVar,
                                         plants$PlantSubSpecies))

## propogate plant key to specimens
specimens$PlantFK <- plants$PlantPK[match(specimens.plant.sp,
                                          plants.plant.sp)]
specimens$UniqueID[is.na(specimens$PlantFK)]

write.csv(dbReadTable(con, 'tblPlant'),
          file='tables/plant.csv', row.names=FALSE)


## *******************************************************
## 5. Pan info:
## *******************************************************

keep <- c('PanColor',
          'PanLocation')
pans <- unique(specimens[keep])
pans$PanPK <- seq_len(nrow(pans))
rownames(pans) <- NULL
dbWriteTable(con, 'tblPan', pans, row.names=FALSE)

## pan key to the specimens table.
pans.pan.info <- paste(pans$PanColor,
                       pans$PanLocation, sep=';')
specimens.pan.info <- paste(specimens$PanColor,
                            specimens$PanLocation, sep=';')
specimens$PanFK <- pans$PanPK[match(specimens.pan.info, pans.pan.info)]

write.csv(dbReadTable(con, 'tblPan'),
          file='tables/pan.csv', row.names=FALSE)

## *******************************************************
## 7. Specimens:
## *******************************************************


keep <- c('UniqueID',
          'TempID',
          'EventID',
          'Collector',
          'SpecimenCount',
          'InsectFK',
          'PlantFK',
          'PanFK',
          'ConditionsFK',
          'GeographyFK')

specimens <- specimens[keep]
specimens <- unique(specimens)
rownames(specimens) <- NULL

dbWriteTable(con, 'tblSpecimens', specimens, row.names=FALSE)

write.csv(dbReadTable(con, 'tblSpecimens'),
          file='tables/specimens.csv', row.names=FALSE)

print(paste("before traditional, dim=", nrow(specimens)))

## *******************************************************
## close connection to database
## *******************************************************
dbDisconnect(con)

## *******************************************************
## make_traditional.R
## *******************************************************

setwd("~/Dropbox/SkyIslands/data/relational/data/relational")
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
drop <- c('GeographyFK')
res.complete <- res.complete[-match(drop, names(res.complete))]

print(paste("after traditional, dim=", nrow(res.complete)))

## set NA to blanks
res.complete[is.na(res.complete)] <- ''

## write full table to csv
write.csv(res.complete, file='traditional/specimens-labels.csv',
          row.names=FALSE)

## **************************************************
## close connection to database
dbDisconnect(con)
