## *******************************************************
## make files for use in relational database
## *******************************************************
rm(list=ls())
setwd("~/Dropbox/skyIslands_saved/data")

## load data
spec.data <- read.csv("raw/specimens.csv",  stringsAsFactors=FALSE,
                      colClasses=c("SampleRound"="character"))

## raw data was entered eith different names than in 2012 but they are
## what we want for the relational database

## *******************************************************
## next we will create the data structures that that we
## will use to construct our relational database
## *******************************************************

spec.data <- spec.data[rep(seq_len(nrow(spec.data)), spec.data$SpecimenCount),]

spec.data$SpecimenCount <- NULL
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


## manually add BBSL numbers - should only need to do this for 2017
## data - will add at point of data entry for future specimens

bbsl <- read.csv("raw/BBSLSpecimens.csv", stringsAsFactors=FALSE)


## only time this will happen is when we get the labels form Terry
spec.data$UniqueID[spec.data$TempID == spec.data$UniqueID] <- bbsl$BarcodeID

write.csv(spec.data, file="relational/original/specimens.csv",
          row.names=FALSE)

source('../../skyIslands/dataPrep/speciesIDs/AssignSpecies.R')


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
write.csv(unique(data.weather), file="relational/original/weather.csv",
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

##Note: when comparing data.geo to H's geography csv:
##  data.geo -> H's Datasheet
##  Site -> NA
##  SubSite -> NA
##  SiteSubSite -> Location_Name
##  MtRange -> NA
##  Forest -> NA (included in Location_Desc)
##  Meadow -> NA (included in Location_Desc)
##  County -> NA
##  State -> NA
##  Country -> NA
##  Lat -> DecimalLat
##  Long -> DecimalLon
##  Elev -> Elev0 (this is the one without m in the cell)

## write unique data to a table
write.csv(data.geo, file="relational/original/geography.csv",
                   row.names=FALSE)
