## *******************************************************
## make files for use in relational database
## *******************************************************
source('dataPrep/src/misc.R')
source('dataPrep/src/taxonomy_gbif.R')
source('dataPrep/relational/src/validate_2017_barcode_map.R')

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

spec.data$Method <- "Net"
spec.data$Method[is.na(spec.data$FieldPlantID) |
                 spec.data$FieldPlantID == ""] <- "Pan"
spec.data$Method[spec.data$SampleRound == 0] <- "Vane"

plant.keys <- read.csv("raw/plants.csv", stringsAsFactors=FALSE)

spec.data$FinalPlantSp <-
    plant.keys$FinalPlantSp[match(spec.data$FieldPlantID,
                                  plant.keys$FieldPlantID)]


## Report unknown field plant codes before they become blanks. We do not
## invent plant identities for those visit specimens.
unknown.plant <- !is.na(spec.data$FieldPlantID) &
  nzchar(trimws(as.character(spec.data$FieldPlantID))) &
  is.na(spec.data$FinalPlantSp)
dir.create('checks/taxonomy', recursive=TRUE, showWarnings=FALSE)
unknown.ids <- sort(unique(as.character(spec.data$FieldPlantID[unknown.plant])))
write.csv(data.frame(FieldPlantID=unknown.ids,
                     NVisitSpecimens=vapply(unknown.ids, function(k)
                       sum(spec.data$FieldPlantID == k, na.rm=TRUE), integer(1))),
    'checks/taxonomy/unmatched_field_plant_codes.csv', row.names=FALSE)
if (any(unknown.plant)) warning(sum(unknown.plant),
  ' specimen visits have unmatched FieldPlantID; see checks/taxonomy.')
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

## Restore historical 2017 barcode numbers only after verifying that the
## original specimen rows, their expansion counts and collecting events have
## NOT shifted. FieldPlantID is an editable field identification; its changes
## are logged below, not treated as changes in specimen order.
spec.data$OriginalUniqueID <- spec.data$UniqueID
spec.data$TempID <- seq_len(nrow(spec.data))
bar.file <- 'raw/2017_legacy_barcode_crosswalk.csv'
if (!file.exists(bar.file)) stop('Missing original 2017 barcode reference: ', bar.file)
bar.map <- read.csv(bar.file, stringsAsFactors=FALSE)
barcode.check <- validate_2017_barcode_map(
  spec.data, bar.map, bbsl,
  audit_file='checks/taxonomy/2017_barcode_field_plant_edits.csv')
rows.2017 <- barcode.check$rows.2017
first.2017 <- barcode.check$first.2017
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

## *******************************************************
## Shared GBIF taxonomy reconciliation: run only AFTER IDs have been
## assigned by AssignSpecies.R. This is the ONLY place where scientific
## names are corrected. Do not alter files under raw/ or speciesIDs/IDs/.
## *******************************************************
check.dir <- 'checks/taxonomy'
dir.create(check.dir, recursive=TRUE, showWarnings=FALSE)
cache.file <- file.path(check.dir, 'gbif_backbone_cache.csv')
## Version-controlled, approved spelling corrections are applied in every run.
## Additional case-by-case approvals can be added in checks/taxonomy/.
overrides.file <- c(
    'dataPrep/speciesIDs/taxon_overrides_approved.csv',
    file.path(check.dir, 'taxon_overrides_approved.csv'))
## For a deliberate full GBIF recheck: set refresh=TRUE below.

identified <- read.csv('relational/original/specimens.csv',
                       stringsAsFactors=FALSE)
identified$InsectOriginalName <- clean_taxon_string(paste(
  ifelse(is.na(identified$Genus), '', identified$Genus),
  ifelse(is.na(identified$Species), '', identified$Species),
  ifelse(is.na(identified$SubSpecies), '', identified$SubSpecies)))

## Collect every distinct name across all four plant inputs, even when
## it appears in only one of them. The same lookup is applied everywhere.
insect.names <- identified$InsectOriginalName
plant.names <- c(plant.keys$FinalPlantSp, identified$FinalPlantSp,
                 veg.data$PlantGenusSpecies, bloom.data$PlantGenusSpecies)
all.names <- c(insect.names, plant.names)
all.kingdoms <- c(rep('Animalia', length(insect.names)),
                  rep('Plantae', length(plant.names)))
all.families <- c(identified$Family,
                  rep('', length(plant.names)))
lookup <- reconcile_taxon_names(all.names, all.kingdoms, all.families,
    cache_file=cache.file, overrides_file=overrides.file, refresh=FALSE)

## Keep an actionable, stable review table of all input name mappings.
taxon.counts <- table(paste(all.kingdoms, clean_taxon_string(all.names),
                            sep='\r'))
lookup$NSourceRecords <- as.integer(taxon.counts[
  paste(lookup$Kingdom, lookup$OriginalName, sep='\r')])
lookup$NSourceRecords[is.na(lookup$NSourceRecords)] <- 0L
lookup$VisitNetworkName <- ifelse(lookup$Kingdom == 'Plantae',
    si_plant_lump(lookup$CorrectedName), lookup$CorrectedName)
write.csv(lookup, file.path(check.dir, 'taxon_name_mapping.csv'),
          row.names=FALSE, na='')
write.csv(lookup[lookup$Status %in%
    c('REVIEW_GBIF_MATCH','REVIEW_FAMILY_CONFLICT','UNCHECKED'), ,drop=FALSE],
    file.path(check.dir, 'taxa_needing_review.csv'),
    row.names=FALSE, na='')
write.csv(lookup[lookup$Status == 'PROVISIONAL_PRESERVED', ,drop=FALSE],
    file.path(check.dir, 'morphospecies_preserved.csv'),
    row.names=FALSE, na='')

## Correct insects once, before tblInsect is constructed. Preserve the
## original label for traceability; leave provisional names untouched.
corrected <- replace_from_taxonomy(identified$InsectOriginalName,
                                  'Animalia', lookup)
identified$InsectGBIFStatus <- lookup$Status[match(
  paste('Animalia', identified$InsectOriginalName, sep='\r'),
  paste(lookup$Kingdom, lookup$OriginalName, sep='\r'))]
identified$InsectGBIFStatus[is.na(identified$InsectGBIFStatus)] <-
  'NO_INSECT_ID'
changed <- nzchar(corrected) & corrected != identified$InsectOriginalName
parts <- strsplit(corrected[changed], ' ', fixed=TRUE)
identified$Genus[changed] <- vapply(parts, `[`, character(1), 1)
identified$Species[changed] <- vapply(parts, function(z)
  if (length(z) >= 2L) z[2] else '', character(1))
identified$SubSpecies[changed] <- vapply(parts, function(z)
  if (length(z) >= 3L) paste(z[-c(1,2)], collapse=' ') else '', character(1))

## Correct plant name references in every input; no downstream lookup
## of normalized_plants.csv and no raw file edits.
plant.keys$FinalPlantSp <- replace_from_taxonomy(
  plant.keys$FinalPlantSp, 'Plantae', lookup)
## Existing analysis convention: field-visit Erigeron IDs are collapsed
## across years to 'Erigeron spp.' Vegetation and bloom observations retain
## species-level identities so floral species diversity is not reduced.
identified$FinalPlantSp <- si_plant_lump(replace_from_taxonomy(
  identified$FinalPlantSp, 'Plantae', lookup))
veg.data$PlantGenusSpecies <- replace_from_taxonomy(
  veg.data$PlantGenusSpecies, 'Plantae', lookup)
bloom.data$PlantGenusSpecies <- replace_from_taxonomy(
  bloom.data$PlantGenusSpecies, 'Plantae', lookup)

## tblPlant must contain every study-level plant taxon used by the
## specimen/veg/bloom source tables, not only raw/plants.csv values.
plant.all <- sort(unique(c(plant.keys$FinalPlantSp,
  identified$FinalPlantSp, veg.data$PlantGenusSpecies,
  bloom.data$PlantGenusSpecies)))
plant.all <- plant.all[nzchar(plant.all)]
plant.lut <- lookup[lookup$Kingdom == 'Plantae', , drop=FALSE]
family.map <- plant.lut$GBIFFamily[match(plant.all, plant.lut$CorrectedName)]
family.map[is.na(family.map)] <- ''
plant.table <- data.frame(PlantGenusSpecies=plant.all,
                          PlantFamily=family.map, stringsAsFactors=FALSE)
write.csv(plant.table, 'relational/original/plant_taxa.csv', row.names=FALSE)
write.csv(identified, 'relational/original/specimens.csv', row.names=FALSE)
write.csv(veg.data, 'relational/original/veg.csv', row.names=FALSE)
write.csv(bloom.data, 'relational/original/bloom.csv', row.names=FALSE)
message('GBIF taxonomy reconciled in 1prep: ', sum(changed),
        ' insect specimen names adjusted; ', length(plant.all),
        ' plant taxon labels represented; see ', check.dir)

