## *******************************************************
## make files for use in relational database
## *******************************************************
rm(list=ls())
setwd("~/Dropbox/SkyIslands/data")


## load data
D <- read.csv("specimenData/cleaned/specimens.csv")

## *******************************************************
## next we will create the data structures that that we
## will use to construct our relational database
## *******************************************************
## create insect specimen database
data <- data.frame(UniqueID=D$temp.id,
                   TempID=D$number,
                   GenID=D$genID,
                   Order=D$order,
                   Family=D$family,
                   Genus=D$genus,
                   SubGenus=D$subgenus,
                   Species=D$species,
                   SubSpecies=D$subspecies,
                   Sex =D$sex,
                   Site=D$site,
                   Collector=D$collector,
                   SampleDay=D$samplingDay,
                   SampleRound=D$samplingRound,
                   Date=D$date,
                   Dos=D$dos,
                   FieldPlantID=D$fieldPlantID,
                   FinalPlantSp=D$finalPlantID,
                   PlantFamily=D$plantFamily,
                   DateDetermined=D$dateDetermined,
                   Determiner=D$determiner,
                   Author=D$author
                   )

check.data <- aggregate(data$Date, list(site = data$Site, day =
                                           data$SampleRound, date =
                                           data$Date), length)

write.csv(data, file="relational/data/original/specimens.csv",
                   row.names=FALSE)
## *******************************************************

## *******************************************************
## create conditions file

W <- read.csv("WeatherConditions.csv")
data <- data.frame(Site=W$site,
                   SampleRound=W$samplingDay,
                   Date=W$date,
                   Dos=W$dos,
                   StartTime=W$startTime,
                   EndTime=W$endTime,
                   TempStart=W$startTemp,
                   TempEnd=W$endTemp,
                   WindStart=W$startWind,
                   WindEnd=W$endWind,
                   SkyStart=W$startWeather,
                   SkyEnd=W$endWeather)

check.data <- aggregate(data$StartTime, list(site = data$Site, day =
                                           data$SampleRound, date =
                                           data$Date), length)

## write unique data to a table
write.csv(unique(data), file="relational/data/original/weather.csv",
                   row.names=FALSE)
## *******************************************************


G <- read.csv("Geography/Geography.csv")
data <- data.frame(Site=G$site,
                   Meadow=G$meadow,
                   MtRange=G$MtRange,
                   Forest=G$forest,
                   County=G$county,
                   State=G$state,
                   Country=G$country,
                   Lat=G$lat,
                   Long=G$long
                   )


## write unique data to a table
write.csv(unique(data), file="relational/data/original/geography.csv",
                   row.names=FALSE)
