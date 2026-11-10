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

##  Preserve provenance before expanding historical
## specimen counts; TempID numbering remains the same as in the original code.
if (any(is.na(spec.data$SpecimenCount) |
        spec.data$SpecimenCount < 0 |
        spec.data$SpecimenCount != floor(spec.data$SpecimenCount)))
  stop('SpecimenCount must contain nonnegative whole numbers.')
if (!'rawRowNum' %in% names(spec.data))
  spec.data$rawRowNum <- seq_len(nrow(spec.data))
spec.data$RawSpecimenCount <- spec.data$SpecimenCount
spec.data <- spec.data[rep(seq_len(nrow(spec.data)),
                           spec.data$SpecimenCount),]
spec.data$ExpandedFromCount <- spec.data$RawSpecimenCount > 1
spec.data$SpecimenCount <- NULL

spec.data$SiteSubSite <-  paste0(spec.data$Site, spec.data$SubSite)

print(paste("original number of specimens", nrow(spec.data)))

## BEWARE CHECK DATE FORMAT
spec.date.format <- "%m/%d/%y"
print(spec.date.format)
print(spec.data$Date[1])
spec.data$Date <- as.Date(spec.data$Date, spec.date.format)

## add subsite one to all 2012 data
spec.2012 <- spec.data$Date < "2013-01-01"
spec.after.2012 <- spec.data$Date > "2013-01-01"

spec.data$Method <- "Net"
spec.data$Method[is.na(spec.data$FieldPlantID) |
                 spec.data$FieldPlantID == ""] <- "Pan"
spec.data$Method[spec.data$SampleRound == 0] <- "Vane"

plant.keys <- read.csv("raw/plants.csv", stringsAsFactors=FALSE)

spec.data$FinalPlantSp <-
    plant.keys$FinalPlantSp[match(spec.data$FieldPlantID,
                                  plant.keys$FieldPlantID)]


spec.data$FinalPlantSp[is.na(spec.data$FinalPlantSp)] <- ""


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

## Use the frozen original order to guard against a new BBSL export with
## inserted rows changing every later specimen ID.
spec.data$OriginalUniqueID <- spec.data$UniqueID
spec.data$TempID <- seq_len(nrow(spec.data))
bar.file <- '../../skyIslands_saved/data/raw/2017_legacy_barcode_crosswalk.csv'
if (!file.exists(bar.file)) stop('Missing original 2017 barcode reference: ', bar.file)
bar.map <- read.csv(bar.file, stringsAsFactors=FALSE)
rows.2017 <- which(spec.data$Year == '2017')
if (length(rows.2017) != 20359L || nrow(bar.map) != 9330L)
  stop('2017 specimen count or historical barcode count has changed. Review IDs first.')
if (anyDuplicated(bar.map$TempID) || anyDuplicated(bar.map$BarcodeID) ||
    anyNA(bar.map$TempID) || anyNA(bar.map$BarcodeID))
  stop('The 2017 barcode reference contains duplicate or missing IDs.')
first.2017 <- rows.2017[seq_len(nrow(bar.map))]
if (!identical(as.integer(spec.data$TempID[first.2017]),
               as.integer(bar.map$TempID)))
  stop('Historical 2017 TempID positions changed; do not shift barcodes.')
## Validate the raw-record provenance associated with
## each historical ID, but do not require collection metadata in BBSL export.
source.fields <- c(SourceRawRowNum='rawRowNum', SourceDate='Date',
  SourceSite='Site', SourceSubSite='SubSite',
  SourceCollector='Collector', SourceFieldPlantID='FieldPlantID',
  SourceNetNumber='NetNumber', SourceSampleRound='SampleRound',
  SourceSpecimenCount='RawSpecimenCount')
normalise.id.field <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c('', 'NA')] <- ''
  x
}
for (field in names(source.fields)) {
  if (!field %in% names(bar.map) ||
      !source.fields[[field]] %in% names(spec.data))
    stop('Missing ID audit field: ', field)
  expected <- normalise.id.field(bar.map[[field]])
  observed <- normalise.id.field(
    spec.data[[source.fields[[field]]]][first.2017])
  if (!identical(expected, observed))
    stop('2017 specimen order differs at ', field,
         '; barcode mapping intentionally not applied.')
}
if (!identical(as.character(bbsl$BarcodeID[seq_len(nrow(bar.map))]),
               as.character(bar.map$BarcodeID)))
  stop('BBSL inventory is not in historical order. Keep raw/BBSLSpecimens.csv unchanged.')
## Assign only historically labeled 2017 records; preserve all other IDs.
spec.data$UniqueID[rows.2017] <- NA_character_
spec.data$UniqueID[first.2017] <- bar.map$BarcodeID
spec.data$UniqueID[is.na(spec.data$UniqueID) |
                   spec.data$UniqueID == ''] <- ''
no.id <- spec.data$UniqueID == ''
spec.data$UniqueID[no.id] <- as.character(spec.data$TempID[no.id])
if (anyDuplicated(spec.data$UniqueID))
  stop('Duplicate specimen UniqueID after restoring historical numbering.')
##  This ID assignment intentionally retains the original TempID
## ordering. Taxonomic ID scripts and parasite files are left unchanged.

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
weather.data$StartDate <- as.Date(weather.data$StartDate, w.date.format)

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

weather.data$Dos <- NULL

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
                               "Pinaleno", "White"))

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
                       Area=geo$Area_m2)

## write unique data to a table
write.csv(data.geo, file="relational/original/geography.csv",
                   row.names=FALSE)


## *******************************************************
## parasite data prep
## *******************************************************

screenings <- c("Apidae", "AspergillusSpp", "AscosphaeraSpp",
                "ApicystisSpp", "CrithidiaExpoeki", "CrithidiaBombi",
                "CrithidiaSpp", "CrithidiaMellificae",
                "NosemaBombi", "NosemaCeranae")

## 2017-2021 data
source('parasite_postivies/all_parasites.R', chdir = TRUE)

para.data <- data.frame(UniqueID=Apidae,
                        TempID=NA,
                        Apidae=1,
                        AspergillusSpp=NA,
                        AscosphaeraSpp=0,
                        ApicystisSpp=0,
                        CrithidiaExpoeki=0,
                        CrithidiaBombi=0,
                        CrithidiaMellificae=0,
                        CrithidiaSpp=0,
                        NosemaCeranae=0,
                        NosemaBombi=0)

para.data$AscosphaeraSpp[para.data$UniqueID %in%
                              AscosphaeraSpp] <- 1
para.data$ApicystisSpp[para.data$UniqueID %in%
                              ApicystisSpp] <- 1
para.data$CrithidiaExpoeki[para.data$UniqueID %in%
                              CrithidiaExpoeki] <- 1
para.data$CrithidiaBombi[para.data$UniqueID %in%
                              CrithidiaBombi] <- 1
para.data$CrithidiaMellificae[para.data$UniqueID %in%
                              CrithidiaMellificae] <- 1
para.data$CrithidiaSpp[para.data$UniqueID %in%
                              CrithidiaSpp] <- 1
para.data$NosemaCeranae[para.data$UniqueID %in%
                              NosemaCeranae] <- 1

## No specimens tested positive for Nosema bombi.
## Preserve the original zero-positive result for successfully extracted
## individuals (Apidae == 1). Do not infer positives from an absent list;
## individuals without successful extraction remain NA after the check below.
para.data$NosemaBombi[para.data$Apidae == 1] <- 0

## contaminated samples
para.data$ApicystisSpp[para.data$UniqueID %in%
                              ApicystisSppNA] <- NA
para.data$CrithidiaExpoeki[para.data$UniqueID %in%
                              CrithidiaExpoekiNA] <- NA
para.data$CrithidiaBombi[para.data$UniqueID %in%
                              CrithidiaBombiNA] <- NA
para.data$CrithidiaSpp[para.data$UniqueID %in%
                              CrithidiaSppNA] <- NA
para.data$CrithidiaMellificae[para.data$UniqueID %in%
                              CrithidiaSppNA] <- NA


## para.data <- rbind(para.data, parasites)

controls <- para.data[grepl("DNActrl", para.data$UniqueID),]
print("Are DNA extraction controls all zero?")
print(all(apply(controls[,screenings], 1, sum, na.rm=TRUE) == 0))

## remove controls
para.data <- para.data[!para.data$UniqueID %in% controls$UniqueID,]

print("before removing species without Apidae positivies")
print(nrow(para.data))

para.data$TempID <- NULL
##  Apidae=1 means DNA extraction succeeded. Only those
## extractions can have pathogen zeros; unknown/failed extractions are NA.
pathogen.cols <- setdiff(names(para.data), c('UniqueID', 'Apidae'))
failed <- is.na(para.data$Apidae) | para.data$Apidae != 1
para.data[failed, pathogen.cols] <- NA

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
