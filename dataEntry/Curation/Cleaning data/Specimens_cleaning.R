

spec <- read.csv("/Users/laurenponisio/Documents/Sky Islands/Data/Specimens091712.csv")

levels(spec$site)

unique(paste(spec$site, spec$date))

levels(spec$finalPlantID)


unique(paste(spec$fieldPlantID, spec$finalPlantID))

spec$finalPlantID <- as.character(spec$finalPlantID)

spec$finalPlantID[spec$finalPlantID ==  "Ceratium nutans"] <- "Cerastium nutans"
spec$finalPlantID[spec$finalPlantID ==  "Germanuim caespitosum"] <- "Geranium caespitosum"
spec$finalPlantID[spec$finalPlantID ==  "Monarda citriodora austromontana"] <- "Monarda citriodora spp. Austromontana" 
spec$finalPlantID[spec$finalPlantID ==  "Pennelia longifolia"] <- "Pennellia longifolia"
spec$finalPlantID[spec$finalPlantID ==  "Penstemon pseudo parvus" ] <- "Penstemon pseudoparvus" 

spec$finalPlantID[spec$finalPlantID == "Cymopterus lemmonii" & spec$fieldPlantID == "Mountain spray"] <- "Holodiscus dumosus"

spec$finalPlantID[spec$finalPlantID == "Monardella odoratissima" & spec$fieldPlantID == "Mint"] <- "Phacelia heterophylla"

spec$finalPlantID[spec$finalPlantID == "Polygonum bistortoides" & spec$fieldPlantID == "yellow pond flower"] <- "Hypericum scouleri"

spec$finalPlantID[spec$finalPlantID == "Geranium viscosissimum" & spec$fieldPlantID == "gentian"] <- "Gentianopsis detonsa"

spec$finalPlantID <- as.factor(spec$finalPlantID)

levels(spec$genSp)

unique(spec$determiner)

levels(spec$dateDetermined)

write.csv(spec, file="/Users/laurenponisio/Documents/Sky Islands/Data/Specimens091812.csv")