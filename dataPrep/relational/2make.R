## *******************************************************
## create relational database
## *******************************************************
source('dataPrep/src/misc.R')
setwd("../skyIslands_saved/data/relational/relational")

climate <- read.csv("../original/climate.csv", as.is=TRUE)
conditions <- read.csv("../original/weather.csv", as.is=TRUE)
specimens <- read.csv("../original/specimens.csv", as.is=TRUE)
geo <- read.csv("../original/geography.csv", as.is=TRUE)
parasites <- read.csv("../original/parasites.csv", as.is=TRUE)
veg <- read.csv("../original/veg.csv", as.is=TRUE)
bloom <- read.csv("../original/bloom.csv", as.is=TRUE)
## Prepared in 1prep: one shared GBIF-corrected plant vocabulary from
## raw/plants + visits + vegetation + bloom surveys.
plants <- read.csv("../original/plant_taxa.csv", as.is=TRUE)

## check that if there is already a database, remove it
if(file.exists("si.db")) file.remove("si.db")
con <- dbConnect(dbDriver("SQLite"), dbname='si.db')

## *******************************************************
## 1. Geographic infomation
## *******************************************************

keep <- c("Site", "Country", "State", "County",
          "Meadow", "Forest",
          "MtRange", "Lat", "Long", "Elev", "Area")

## Sites are representative site coordinates, NOT individual net locations.
geography <- geo[!duplicated(geo$Site), keep, drop=FALSE]
geography <- geography[order(geography$Site), , drop=FALSE]
if (anyDuplicated(geography$Site) || anyNA(geography$Site))
  stop('Geography must contain one nonmissing record per Site.')

## generate primary geography key
geography <- cbind(GeographyPK=seq_len(nrow(geography)), geography)
rownames(geography) <- NULL
dbWriteTable(con, "tblGeography", geography, row.names=FALSE)

## Propagate geography key to the conditions table.
conditions$GeographyFK <-
    geography$GeographyPK[match(conditions$Site,
                                geography$Site)]

print(paste("conditions without site matches",
            unique(conditions$Site[is.na(conditions$GeographyFK)])))

## Propagate geography key to the specimens table.
specimens$GeographyFK <- geography$GeographyPK[match(specimens$Site,
                                                     geography$Site)]

## print("specimens without site matches",
##             unique(specimens$UniqueID[is.na(specimens$GeographyFK)]))

## Propagate geography key to the veg table.
veg$GeographyFK <- geography$GeographyPK[match(veg$Site,
                                               geography$Site)]

print(paste("veg quads without site matches",
            unique(veg$Site[is.na(veg$GeographyFK)])))


## Propagate geography key to the bloom table.
bloom$GeographyFK <- geography$GeographyPK[match(bloom$Site,
                                               geography$Site)]

print(paste("bloom data without site matches",
            unique(bloom$Site[is.na(bloom$GeographyFK)])))



## Propagate geography key to the climate table.
climate$GeographyFK <- geography$GeographyPK[match(climate$Site,
                                               geography$Site)]

print(paste("climate data without site matches",
            unique(climate$Site[is.na(climate$GeographyFK)])))

## write a .csv version of this table (just for ease of viewing)
write.csv(dbReadTable(con, "tblGeography"),
          file="tables/geography.csv", row.names=FALSE)

dbListTables(con)


print("conditions with no site key")
print(sum(is.na(conditions$GeographyFK)))
no.geo.cond <- conditions[is.na(conditions$GeographyFK),]
write.csv(no.geo.cond, file="../../cleaning/conditions_no_geo_key.csv")

print("specimens with no site key")
print(sum(is.na(specimens$GeographyFK)))
no.geo.spec <- specimens[is.na(specimens$GeographyFK),]
print(no.geo.spec$SampleID)
write.csv(no.geo.spec, file="../../cleaning/specimens_no_geo_key.csv")
unique(no.geo.spec$Stand)

## *******************************************************
## 2. Conditions
## *******************************************************
## 2.1 for insect specimens
## *******************************************************

## CHECK DATE FORMAT
conditions$Date <- as.Date(conditions$Date, format = "%Y-%m-%d")
specimens$Date <- as.Date(specimens$Date, format = "%Y-%m-%d")
veg$Date <- as.Date(veg$Date, format = "%Y-%m-%d")
bloom$Date <- as.Date(bloom$Date, format = "%Y-%m-%d")

conditions$Date  <- as.character(conditions$Date)
specimens$Date  <- as.character(specimens$Date)
veg$Date  <- as.character(veg$Date)
bloom$Date  <- as.character(bloom$Date)

## Remove columns that will become confusing down the line
climate$StartDate <- NULL
climate$EndDate <- NULL

##  Site x Year x SampleRound must be a single
## climate record before joining, or weather rows can multiply. Retain
## the FIRST climate row as the previous merge intended; write conflicts
## for review rather than silently expanding the conditions table.
for (d in c('conditions', 'climate', 'specimens', 'veg', 'bloom')) {
  x <- get(d)
  x$Year <- as.character(x$Year)
  x$SampleRound <- as.character(x$SampleRound)
  assign(d, x)
}
clim.keys <- c('Site', 'Year', 'SampleRound')
clim.duplicated <- duplicated(climate[clim.keys]) |
                   duplicated(climate[clim.keys], fromLast=TRUE)
dir.create('../../checks/relational_audit', recursive=TRUE,
           showWarnings=FALSE)
write.csv(climate[clim.duplicated, , drop=FALSE],
  file='../../checks/relational_audit/climate_duplicate_rounds.csv',
  row.names=FALSE)
if (any(clim.duplicated))
  warning('Repeated climate Site/Year/Round keys: ',
          sum(clim.duplicated), ' rows; first record retained.')
climate <- climate[!duplicated(climate[clim.keys]), , drop=FALSE]
conditions.n <- nrow(conditions)
conditions <- dplyr::left_join(
  conditions,
  dplyr::select(climate, -GeographyFK),
  by = clim.keys
)
if (nrow(conditions) != conditions.n)
  stop('Climate join unexpectedly changed the weather row count.')

## Temporarily identify unique combinations:
keep.spec <- c("GeographyFK", "Date", "Method", "NetNumber")
conditions$cond.code <- apply(conditions[keep.spec], 1, paste, collapse=";")
specimens$cond.code <- apply(specimens[keep.spec], 1, paste, collapse=";")

## make table
## Select variables for the conditions table
keep <- c(
  "Date", "SampleRound", "NetNumber", "Method",
  "StartTime", "EndTime",
  "TempStart", "TempEnd",
  "WindStart", "WindEnd",
  "SkyStart", "SkyEnd",
  "MeanTmean", "MeanTmax", "MeanTmin",
  "CumPrecip", "MeanVPD",
  "Tmean_anom", "Tmax_anom", "Tmin_anom",
  "Precip_anom", "VPD_anom",
  "APi", "GDD",
  "GeographyFK", "cond.code"
)

cond <- unique(conditions[keep])
##  Multiple weather rows can share the same SITE-level
## event key (e.g., measurements from different subsites). The original
## match() assigned the first record. Preserve that behavior, but keep
## exactly ONE dimension row per event and audit conflicting records.
cond.duplicate <- duplicated(cond$cond.code) |
                  duplicated(cond$cond.code, fromLast=TRUE)
write.csv(cond[cond.duplicate, , drop=FALSE],
  file='../../checks/relational_audit/weather_repeated_site_events.csv',
  row.names=FALSE)
if (any(cond.duplicate))
  message('## REVIEW: ', length(unique(cond$cond.code[cond.duplicate])),
          ' site-level weather event keys have multiple measurements;',
          ' retaining the first entry.')
cond <- cond[!duplicated(cond$cond.code), , drop=FALSE]
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

bad <- (specimens[, c("Site", "Date",
                            "NetNumber", "SubSite")][is.na(specimens$ConditionsFK),])

write.csv(bad, file="weather_problems.csv", row.names=FALSE)


bad.spec <- (specimens[, c("UniqueID", "Site", "Date",
                            "NetNumber", "SubSite")][is.na(specimens$ConditionsFK),])

write.csv(bad.spec, file="../../cleaning/spec_weather_problems.csv", row.names=FALSE)

print("specimens without condition keys")
spec.no.con.key <- unique(specimens$UniqueID[is.na(specimens$ConditionsFK)])
print(spec.no.con.key)

write.csv(specimens[is.na(specimens$ConditionsFK),],
          file="../../cleaning/specimens_no_cond_key.csv")

write.csv(unique(specimens[is.na(specimens$ConditionsFK),
                           c("Site", keep.spec)]),
          file="../../cleaning/cond_specimens_no_cond_key.csv")


write.csv(dbReadTable(con, "tblConditions"),
          file="tables/conditions.csv", row.names=FALSE)

## *******************************************************
## 2.1 sample dates  for veg datasets
## *******************************************************

keep <- c("Site", "Year", "SampleRound")
## One sample FK per Site x Year x SampleRound,
## not per visit date, since a round can span multiple days. 
conditions$samp.code <- gsub(' ', '', apply(conditions[keep], 1,
                                           paste, collapse=';'))
veg$samp.code <- gsub(' ', '', apply(veg[keep], 1, paste, collapse=';'))
bloom$samp.code <- gsub(' ', '', apply(bloom[keep], 1, paste, collapse=';'))
keep <- c('Site', 'Date', 'SampleRound', 'samp.code')
samp <- conditions[!duplicated(conditions$samp.code), keep, drop=FALSE]
rownames(samp) <- NULL
samp <- cbind(SamplePK=seq_len(nrow(samp)), samp)

## Don't upload the cond.code column
dbWriteTable(con, "tblSample", samp[-ncol(samp)], row.names=FALSE)

## Propagate conditions key to the veg table.
veg$SampleFK <-
    samp$SamplePK[match(veg$samp.code, samp$samp.code)]

bloom$SampleFK <-
    samp$SamplePK[match(bloom$samp.code, samp$samp.code)]

print(paste("site, round in quadrat data without condition key",
            unique(paste(veg$Site,
                         veg$SampleRound, veg$Date)[is.na(veg$SampleFK)])))

print(paste("site, round in bloom data without condition key",
            unique(paste(bloom$Site,
                         bloom$Date, bloom$SampleRound)[is.na(bloom$SampleFK)])))


## *******************************************************
## 3. Insect species:
## *******************************************************

keep <- c("Order", "Family", "Genus", "SubGenus", "Species",
          "SubSpecies", "Author")

insects <- specimens[keep]
insects <- unique(insects)

insects$gen.sp <- paste(insects$Order,
                        insects$Family,
                        insects$Genus,
                        insects$SubGenus,
                        insects$Species,
                        insects$SubSpecies,
                        insects$Author,
                        ## insects$Sex,
                        ## insects$Determiner,
                        sep=";")
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
                          specimens$SubSpecies,
                          ## UPDATED [CONS-10]: Same fields on both sides.
                          specimens$Author,
                          ## specimens$Sex,
                          ## specimens$Determiner,
                          sep=";")
specimens$InsectFK <- insects$InsectPK[match(specimens$gen.sp,
                                             insects$gen.sp)]

print(paste("insects without insect IDs",
            specimens$UniqueID[is.na(specimens$InsectFK)]))

write.csv(dbReadTable(con, "tblInsect"),
          file="tables/insect.csv", row.names=FALSE)

## *******************************************************
## 4. Plant species:
## *******************************************************

## Preserve a complete scientific name explicitly, instead of
## reassembling it differently for different surveys. Non-binomial labels
## such as "Agoseris spp." and infraspecific names are valid study labels.
plant.names <- trimws(as.character(plants$PlantGenusSpecies))
plant.names <- unique(plant.names[!is.na(plant.names) & nzchar(plant.names)])
if (anyDuplicated(plant.names)) stop('Duplicate canonical plant names.')
parts <- strsplit(plant.names, ' ', fixed=TRUE)
plants <- data.frame(
  PlantPK=seq_along(plant.names),
  PlantGenusSpecies=plant.names,
  PlantGenus=vapply(parts, `[`, character(1), 1),
  PlantSpecies=vapply(parts, function(z) if (length(z) >= 2L) z[2] else '',
                      character(1)),
  ## Keep all tokens after the species epithet together, including
  ## var./subsp. rank markers when present.
  PlantVar=vapply(parts, function(z) if (length(z) > 2L)
                   paste(z[-c(1, 2)], collapse=' ') else '', character(1)),
  PlantSubSpecies='',
  stringsAsFactors=FALSE)
plant.family <- read.csv('../original/plant_taxa.csv',
                         stringsAsFactors=FALSE)
plants$PlantFamily <- plant.family$PlantFamily[
   match(plants$PlantGenusSpecies, plant.family$PlantGenusSpecies)]
rownames(plants) <- NULL
dbWriteTable(con, 'tblPlant', plants, row.names=FALSE)

## Audit names not matched to a plant FK; only truly blank names
## should remain unmatched after collecting all plant sources in 1prep.
plant_fk <- function(x) plants$PlantPK[match(trimws(as.character(x)),
                                             plants$PlantGenusSpecies)]
specimens$PlantFK <- plant_fk(specimens$FinalPlantSp)
veg$PlantFK <- plant_fk(veg$PlantGenusSpecies)
bloom$PlantFK <- plant_fk(bloom$PlantGenusSpecies)
for (nm in c('specimens','veg','bloom')) {
  dat <- get(nm)
  col <- if (nm == 'specimens') 'FinalPlantSp' else 'PlantGenusSpecies'
  unmatched <- is.na(dat$PlantFK) & !is.na(dat[[col]]) &
    nzchar(trimws(as.character(dat[[col]])))
  if (any(unmatched)) stop(nm, ': ', sum(unmatched),
    ' named plant records cannot match the shared plant table: ',
    paste(utils::head(unique(dat[[col]][unmatched]), 10), collapse=', '))
}

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
## 6. Veg quadrat data
## *******************************************************

keep <- c('Quadrat')
quads <- unique(veg[keep])

quads$QuadPK <- seq_len(nrow(quads))

rownames(quads) <- NULL
dbWriteTable(con, 'tblQuad', quads, row.names=FALSE)

write.csv(dbReadTable(con, 'tblQuad'),
          file='tables/quads.csv', row.names=FALSE)

veg$QuadFK <- quads$QuadPK[match(veg$Quadrat, quads$Quadrat)]

## *******************************************************
## 7. Specimens:
## *******************************************************

##  original event fields only as provenance.
## These never override or substitute for a missing ConditionsFK.
specimens$SourceDate <- specimens$Date
specimens$SourceSite <- specimens$Site
specimens$SourceSampleRound <- specimens$SampleRound
specimens$SourceMethod <- specimens$Method
specimens$SourceNetNumber <- specimens$NetNumber

keep <- c('UniqueID',
          'TempID',
    #      'EventID',
          'Collector',
          'InsectOriginalName', 'InsectGBIFStatus',
          'SourceDate', 'SourceSite', 'SourceSampleRound',
          'SourceMethod', 'SourceNetNumber',
          'InsectFK',
          'PlantFK',
          'PanFK',
          'ConditionsFK',
          'GeographyFK',
          'Determiner',
          'Sex')

specimens <- specimens[keep]
specimens <- unique(specimens)
rownames(specimens) <- NULL

## *******************************************************
## 8. Add parasite data:
## *******************************************************
##  Parasite table must contain at most one result
## per UniqueID, or merge() silently multiplies specimens.
if (anyDuplicated(parasites$UniqueID)) {
  write.csv(parasites[duplicated(parasites$UniqueID) |
                      duplicated(parasites$UniqueID, fromLast=TRUE), , drop=FALSE],
    '../../checks/relational_audit/duplicate_parasite_ids.csv',
    row.names=FALSE)
  stop('Duplicate parasite UniqueID; see checks/relational_audit.')
}
specimens.n <- nrow(specimens)
specimens <- merge(specimens, parasites, by='UniqueID',
                   all.x=TRUE, sort=FALSE)
if (nrow(specimens) != specimens.n)
  stop('Parasite join changed the number of specimens.')
if (anyDuplicated(specimens$UniqueID))
  stop('Final specimens contain duplicate UniqueID values.')


## *******************************************************
## 9. write final specimen table and close connection to database
## *******************************************************

dbWriteTable(con, 'tblSpecimens', specimens, row.names=FALSE)

write.csv(dbReadTable(con, 'tblSpecimens'),
          file='tables/specimens.csv', row.names=FALSE)

print(paste("before specimen traditional, dim=", nrow(specimens)))

## *******************************************************
## 10. quadrats:
## *******************************************************

##  Retain actual vegetation survey date.
keep <- c('QuadFK',
          'PlantFK',
          'SampleFK',
          'Date',
          'BloomStatus',
          'PlantCount',
          'NumBlooms',
          'Notes')

veg <- veg[keep]
veg <- unique(veg)

dbWriteTable(con, 'tblVeg', veg, row.names=FALSE)

write.csv(dbReadTable(con, 'tblVeg'),
          file='tables/veg.csv', row.names=FALSE)

print(paste("before veg traditional, dim=", nrow(veg)))

## *******************************************************
## 11. Bloom data:
## *******************************************************

##  Retain actual bloom survey date.
keep <- c('PlantFK',
          'SampleFK',
          'Date',
          'BloomStatus',
          "NumBlooms",
          'Notes')

bloom <- bloom[keep]
bloom <- unique(bloom)

dbWriteTable(con, 'tblBloom', bloom, row.names=FALSE)

write.csv(dbReadTable(con, 'tblBloom'),
          file='tables/bloom.csv', row.names=FALSE)

print(paste("before bloom traditional, dim=", nrow(bloom)))

dbDisconnect(con)

