## *******************************************************
## make files for use in relational database
## *******************************************************
source('dataPrep/src/misc.R')

setwd("../skyIslands_saved/data")
## load data
spec.data <- read.csv("raw/specimens.csv",  stringsAsFactors=FALSE,
                      colClasses=c("SampleRound"="character"))

## raw data was entered with different names than in 2012 but they are
## what we want for the relational database

## *******************************************************
## next we will create the data structures that that we
## will use to construct our relational database
## *******************************************************

spec.data <- spec.data[rep(seq_len(nrow(spec.data)),
                           spec.data$SpecimenCount),]

spec.data$SpecimenCount <- NULL

spec.data$SiteSubSite <-  paste0(spec.data$Site, spec.data$SubSite)
## spec.data$Site <- gsub("[0-9]", "", spec.data$Site)

print(paste("original number of specimens", nrow(spec.data)))

## BEWARE CHECK DATE FORMAT
spec.date.format <- "%m/%d/%y"
print(spec.date.format)
print(spec.data$Date[1])
spec.data$Date <- as.Date(spec.data$Date, spec.date.format)

## add subsite one to all 2012 data
spec.2012 <- spec.data$Date < "2013-01-01"
spec.after.2012 <- spec.data$Date > "2013-01-01"

## spec.data$SiteSubSite[spec.2012] <-
##     paste0(spec.data$SiteSubSite[spec.2012], "1")

spec.data$Method <- "Net"
spec.data$Method[is.na(spec.data$FieldPlantID) |
                 spec.data$FieldPlantID == ""] <- "Pan"
spec.data$Method[spec.data$SampleRound == 0] <- "Vane"

plant.keys <- read.csv("raw/plants.csv", stringsAsFactors=FALSE)

spec.data$FinalPlantSp <-
    plant.keys$FinalPlantSp[match(spec.data$FieldPlantID,
                                  plant.keys$FieldPlantID)]

## spec.data$UniqueID[spec.2017] <- spec.data$TempID[spec.2017] <-
##     1:length(spec.data$UniqueID[spec.2017])


spec.data$FinalPlantSp[is.na(spec.data$FinalPlantSp)] <- ""

## spec.data$NetNumber <-
##     sapply(strsplit(spec.data$SampleRound, ".", fixed=TRUE),
##            function(x) x[2])

## spec.data$SampleRound <-
##     sapply(strsplit(spec.data$SampleRound, ".", fixed=TRUE),
##            function(x) x[1])


check.data.spec <- aggregate(spec.data$Date,
                             list(site = spec.data$Site,
                                  date =
                                      spec.data$Date,
                                  subsite=spec.data$SubSite,
                                  netnumber=spec.data$NetNumber,
                                  sampleround=spec.data$SampleRound
                                  ),
                             length)

check.data.spec <- check.data.spec[order(check.data.spec$date),]
write.csv(check.data.spec, file="spec_counts.csv", row.names=FALSE)

## manually add BBSL numbers - should only need to do this for 2017
## data - will add at point of data entry for future specimens

bbsl <- read.csv("raw/BBSLSpecimens.csv", stringsAsFactors=FALSE)


## only time this will happen is when we get the labels form Terry
spec.data$UniqueID[spec.data$Year == "2017"] <-
    bbsl$BarcodeID[1:sum(spec.data$Year == "2017")]

spec.data$TempID <- 1:nrow(spec.data)
spec.data$UniqueID[is.na(spec.data$UniqueID)] <-  ""

spec.data$UniqueID[spec.data$UniqueID == ""] <-
    spec.data$TempID[spec.data$UniqueID == ""]

## spec.data$UniqueID[20794:nrow(spec.data)] <- 20794:nrow(spec.data)

write.csv(spec.data, file="relational/original/specimens.csv",
          row.names=FALSE)

source('../../skyIslands/dataPrep/speciesIDs/AssignSpecies.R')


## *******************************************************
## create conditions file
weather.data <- read.csv("raw/weather.csv")

weather.data$SiteSubSite <-  paste0(weather.data$Site,
                                    weather.data$SubSite)

## BEWARE CHECK DATE
w.date.format <- "%m/%d/%y"
print(w.date.format)
print(weather.data$Date[1])
weather.data$Date <- as.Date(weather.data$Date, w.date.format)

convertFtoC <- function(x){
  y <- round((x-32)*(5/9),2)
  return(y)
}

weather.data$TempStart[weather.data$TempStart > 50  & is.finite(weather.data$TempStart)] <- 
  convertFtoC(weather.data$TempStart[weather.data$TempStart > 50 & is.finite(weather.data$TempStart)])


weather.data$TempEnd[weather.data$TempEnd == ""] <- NA
weather.data$TempEnd <- as.numeric(weather.data$TempEnd)
weather.data$TempEnd[weather.data$TempEnd > 50 & is.finite(weather.data$TempEnd)] <- 
  convertFtoC(weather.data$TempEnd[weather.data$TempEnd > 50 & is.finite(weather.data$TempEnd)])


check.weather.data <- aggregate(weather.data$StartTime,
                                list(site = weather.data$Site,
                                     day = weather.data$SampleRound,
                                     date = weather.data$Date), length)

## write unique data to a table
write.csv(unique(weather.data), file="relational/original/weather.csv",
                   row.names=FALSE)
## *******************************************************

geo <- read.csv("raw/geography.csv")

colnames(geo) <- gsub("\\.", "", colnames(geo))

## geo$Site <- gsub("[0-9]", "", geo$Site)

geo$SubSite <- substr(geo$SiteSubSite, 3, 3)

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
                       SiteSubSite=geo$SiteSubSite,
                       MtRange=geo$MtRange,
                       Forest=geo$Forest,
                       Meadow=geo$Meadow,
                       County=geo$County,
                       State=geo$State,
                       Country=geo$Country,
                       Lat=geo$DecimalLat,
                       Long=geo$DecimalLon,
                       Elev=geo$Elev0,
                       Area=geo$Area)

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


## *******************************************************
## parasite data prep
## *******************************************************

screenings <- c("Apidae", "AspergillusSpp", "AscosphaeraSpp",
                "ApicystisSpp", "CrithidiaExpoeki", "CrithidiaBombi",
                "CrithidiaSpp",
                "NosemaBombi", "NosemaCeranae")

## con
## par.cols <-  c("UniqueID",
##               "TempID",
##               "Apidae",
##               "AspergillusSpp",
##               "AscosphaeraSpp",
##               "ApicystisSpp",
##               "CrithidiaExpoeki",
##               "CrithidiaBombi",
##               "CrithidiaSpp",
##               "NosemaCeranae",
##               "NosemaBombi")
## ## 2018 data
## para.data <- read.csv("raw/parasites.csv",  stringsAsFactors=FALSE)
## para.data$CrithidiaSpp <- NA
## para.data[, par.cols)]

## dir.pars <- "parasite_postivies/indiv_parasites"

## for(i in par.cols[-c(1,2)]){
##   positives <- para.data$UniqueID[para.data[,i] == 1]
 
##   write.table(glue_collapse(sort(positives), sep=", "),
##               file=file.path(dir.pars, sprintf("/%s_2018.txt",
##                                                i)), sep=",",
##               row.names=FALSE)
## }


## 2021-2022 data
source('parasite_postivies/all_parasites.R', chdir = TRUE)

para.data <- data.frame(UniqueID=Apidae,
                             TempID=NA,
                             Apidae=1,
                             AspergillusSpp=NA,
                             AscosphaeraSpp=0,
                             ApicystisSpp=0,
                             CrithidiaExpoeki=0,
                             CrithidiaBombi=0,
                             CrithidiaSpp=0,
                             NosemaCeranae=0, ## change to 0 after
                             ## the screenings are completed
                             NosemaBombi=0)

para.data$AscosphaeraSpp[para.data$UniqueID %in%
                              AscosphaeraSpp] <- 1
para.data$ApicystisSpp[para.data$UniqueID %in%
                              ApicystisSpp] <- 1
para.data$CrithidiaExpoeki[para.data$UniqueID %in%
                              CrithidiaExpoeki] <- 1
para.data$CrithidiaBombi[para.data$UniqueID %in%
                              CrithidiaBombi] <- 1
para.data$CrithidiaSpp[para.data$UniqueID %in%
                              CrithidiaSpp] <- 1

## contaminated samples
para.data$ApicystisSpp[para.data$UniqueID %in%
                              ApicystisSppNA] <- NA
para.data$CrithidiaExpoeki[para.data$UniqueID %in%
                              CrithidiaExpoekiNA] <- NA
para.data$CrithidiaBombi[para.data$UniqueID %in%
                              CrithidiaBombiNA] <- NA
para.data$CrithidiaSpp[para.data$UniqueID %in%
                              CrithidiaSppNA] <- NA


## para.data <- rbind(para.data, parasites)

controls <- para.data[grepl("DNActrl", para.data$UniqueID),]
print("Are DNA extraction controls all zero?")
print(all(apply(controls[,screenings], 1, sum, na.rm=TRUE) == 0))

## remove controls
para.data <- para.data[!para.data$UniqueID %in% controls$UniqueID,]

print("before removing species without Apidae positivies")
print(nrow(para.data))

## para.data <- para.data[para.data$Apidae == 1,]
## print("after removing species without Apidae positivies")
## print(nrow(para.data))

para.data$TempID <- NULL
## para.data$Apidae <- NULL

write.csv(para.data, file="relational/original/parasites.csv",
                   row.names=FALSE)

## *******************************************************
## veg
## *******************************************************

veg.data <- read.csv("raw/quadratveg.csv",  stringsAsFactors=FALSE)

## there were a lot of strange quadrat numbers in the data which are
## not possible. We only had 3 subsites, so subsite 4-5 are not
## possible.

veg.data$Quadrat[veg.data$Quadrat == "C4"]  <- "C3"
veg.data$Quadrat[veg.data$Quadrat == "C5"]  <- "C3"

veg.data$Quadrat[veg.data$Quadrat == "E4"]  <- "E3"
veg.data$Quadrat[veg.data$Quadrat == "E5"]  <- "E3"

veg.data$Quadrat[veg.data$Quadrat == "N4"]  <- "N3"
veg.data$Quadrat[veg.data$Quadrat == "N5"]  <- "N3"

veg.data$Quadrat[veg.data$Quadrat == "W4"]  <- "W3"
veg.data$Quadrat[veg.data$Quadrat == "W5"]  <- "W3"

veg.data$Quadrat[veg.data$Quadrat == "S4"]  <- "S3"
veg.data$Quadrat[veg.data$Quadrat == "S5"]  <- "S3"

veg.data$Quadrat[veg.data$Quadrat == "ME4"]  <- "ME3"
veg.data$Quadrat[veg.data$Quadrat == "ME5"]  <- "ME3"

veg.data$Quadrat[veg.data$Quadrat == "MN4"]  <- "ME3"
veg.data$Quadrat[veg.data$Quadrat == "MN5"]  <- "ME3"

veg.data$Quadrat[veg.data$Quadrat == "MS4"]  <- "MS3"
veg.data$Quadrat[veg.data$Quadrat == "MS5"]  <- "MS3"

veg.data$Quadrat[veg.data$Quadrat == "MW4"]  <- "MW3"
veg.data$Quadrat[veg.data$Quadrat == "MW5"]  <- "MW3"


veg.data$Quadrat[veg.data$Quadrat == "CW2"]  <- "MC2"


## cannot have middle center....
veg.data$Quadrat[veg.data$Quadrat == "MC2"]  <- "ME2"
veg.data$Quadrat[veg.data$Quadrat == "NC1"]  <- "C1"


# reversed order of names
veg.data$Quadrat[veg.data$Quadrat == "WM1"]  <- "MW1"
veg.data$Quadrat[veg.data$Quadrat == "WM2"]  <- "MW2"
veg.data$Quadrat[veg.data$Quadrat == "NM3"]  <- "MN3"


veg.data$SubSite <- parse_number(veg.data$Quadrat)

## lots of site name issues
veg.data$Site[veg.data$Site == "MG"]  <- "MM"
veg.data$Site[veg.data$Site == "BP"]  <- "CH"
veg.data$Site[veg.data$Site == "SD"]  <- "SC"

veg.data$SiteSubSite <-  paste0(veg.data$Site,
                                veg.data$SubSite)

veg.data$PlantGenusSpecies <-
    fix.white.space(veg.data$PlantGenusSpecies)

veg.data$PlantGenusSpecies[veg.data$PlantGenusSpecies == "NA"] <- ""

## CHECK DATE FORMAT
veg.data$Date <- as.Date(veg.data$Date, "%m/%d/%y")

write.csv(veg.data, file="relational/original/veg.csv",
                   row.names=FALSE)



## *******************************************************
## blooms
## *******************************************************

bloom.data <- read.csv("raw/blooms.csv",  stringsAsFactors=FALSE)


bloom.data$SiteSubSite <-  paste0(bloom.data$Site,
                                bloom.data$SubSite)

bloom.data$PlantGenusSpecies <-
    fix.white.space(bloom.data$PlantGenusSpecies)

bloom.data$PlantGenusSpecies[bloom.data$PlantGenusSpecies == "NA"] <-
    ""

bloom.data$NumBlooms <- bloom.data$NumFlower
bloom.data$NumFlower <- NULL



## CHECK DATE FORMAT
bloom.data$Date <- as.Date(bloom.data$Date, "%m/%d/%y")

## sort(unique(bloom.data$PlantGenusSpecies))

write.csv(bloom.data, file="relational/original/bloom.csv",
                   row.names=FALSE)
