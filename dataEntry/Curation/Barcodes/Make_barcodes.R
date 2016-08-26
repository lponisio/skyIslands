spec <- read.csv("/Users/laurenponisio/Documents/Sky Islands/Data/BarcodesSpecimens091812.csv")

spec$uniqueID <- paste(spec$essig, spec$project, spec$site, spec$date, spec$number, sep='_')

write.csv(spec, "/Users/laurenponisio/Documents/Sky Islands/Data/BarcodesSpecimens091912_2.csv", row.names= FALSE)

