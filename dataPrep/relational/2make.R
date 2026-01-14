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
plants <- read.csv("../../raw/plants.csv",
                   as.is=TRUE)

## check that if there is already a database, remove it
if(file.exists("si.db")) file.remove("si.db")
con <- dbConnect(dbDriver("SQLite"), dbname='si.db')

## *******************************************************
## 1. Geographic infomation
## *******************************************************

keep <- c("Site", "Country", "State", "County",
          "Meadow", "Forest",
          "MtRange", "Lat", "Long", "Elev", "Area")

geography <- unique(geo[keep])
## next sort into alphabetical order
geography <- geography[match(sort(geography$Site),
                             geography$Site),]

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

## Merge climate to conditions. This is not an ideal implementation
## because the climate data will be repeated for SampleRounds that
## endend multiple days. But there isn't a clean way to merge the
## different dimensions of data without creating another redundant
## key. 
conditions <- merge(conditions, climate, all.x=TRUE)

## Temporarily identify unique combinations:
keep.spec <- c("GeographyFK", "Date", "Method", "NetNumber")
conditions$cond.code <- apply(conditions[keep.spec], 1, paste, collapse=";")
specimens$cond.code <- apply(specimens[keep.spec], 1, paste, collapse=";")

## make table
keep <- c("Date", "SampleRound", "NetNumber", "Method", "StartTime",
          "EndTime", "TempStart", "TempEnd", "WindStart", "WindEnd",
          "SkyStart", "SkyEnd", "SpringPrecip", "CumulativePrecip",
          "RoundPrecip", "SpringTmean", "RoundTmean",
          "CumulativeTmean", "SpringTmeanAnom", "RoundTmeanAnom",
          "CumulativeTmeanAnom",
          "GeographyFK", "cond.code")

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
conditions$samp.code <- gsub(" ", "", apply(conditions[keep], 1, paste, collapse=";"))
veg$samp.code <- apply(veg[keep], 1, paste, collapse=";")
bloom$samp.code <- apply(bloom[keep], 1, paste, collapse=";")

keep <- c("Site", "Date", "SampleRound", "samp.code")
samp <- unique(conditions[keep])
rownames(samp) <- NULL
samp <- cbind(SamplePK=seq_len(nrow(samp)), samp)

## Don't upload the cond.code column
dbWriteTable(con, "tblSample", samp[-ncol(samp)], row.names=FALSE)

## ## Propagate conditions key to the conditions table.
## conditions$SampleFK <-
##     samp$SamplePK[match(conditions$samp.code, samp$samp.code)]

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
                         bloom$Date, veg$Date)[is.na(bloom$SampleFK)])))


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

keep <- c('FinalPlantSp')
plants <- fix.white.space(paste(plants$PlantGenus,
                                         plants$PlantSpecies,
                                         plants$PlantSubspecies))
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

print(paste("spec IDs without plants IDs",
            specimens$UniqueID[is.na(specimens$PlantFK)]))

unique(print(paste("spec without plants IDs",
            unique(specimens$FinalPlantSp[is.na(specimens$PlantFK)]))))

# *********
## Propagate plant key to the veg table.
veg$PlantFK <- plants$PlantPK[match(veg$PlantGenusSpecies,
                                          plants.plant.sp)]

print(paste("veg without plants IDs",
            unique(veg$PlantGenusSpecies[is.na(veg$PlantFK)])))

# *********
## Propagate plant key to the bloom table.
bloom$PlantFK <- plants$PlantPK[match(bloom$PlantGenusSpecies,
                                          plants.plant.sp)]

print(paste("bloom without plants IDs",
            unique(bloom$PlantGenusSpecies[is.na(bloom$PlantFK)])))


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

keep <- c('UniqueID',
          'TempID',
    #      'EventID',
          'Collector',
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
specimens <- merge(specimens, parasites, by="UniqueID",
                   all.x=TRUE)


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

keep <- c('QuadFK',
          'PlantFK',
          'SampleFK',
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

keep <- c('PlantFK',
          'SampleFK',
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

