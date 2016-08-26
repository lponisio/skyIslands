rm(list=ls())
setwd("~/Dropbox/SI_data_entry")
spec <- read.csv("Specimens03122013.csv")

source('Species_ids/Butterfly_morpho (1).R')
geo <- read.csv("Geography.csv")
## modify date:
spec$date <- paste("0", spec$date, sep="")

## create unique ID:
spec$uniqueID <- paste(rep("EMEC", nrow(spec)),
                       rep("SI", nrow(spec)),
                       spec$site,
                       spec$date,
                       spec$number,
                       sep='_')
## add column with species type:
findID <- paste(spec$site, spec$date, spec$number, sep='_')

type <- rep(names(butterfly.list), sapply(butterfly.list, length))
matched <- match(findID, unlist(butterfly.list))
pos <- which(!is.na(matched))
key <- rep("", nrow(spec))
key[pos] <- type[matched[pos]]
spec$key <- key

butterfly.table <- read.csv("species_table.csv")

merge.spec <- merge(spec, butterfly.table, by="key")

merge.spec.butterfly <- merge.spec[!is.na(merge.spec$genus_species),]

write.csv(merge.spec.butterfly, "butterfly_labels_0409.csv")