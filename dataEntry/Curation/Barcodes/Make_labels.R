spec <- read.csv("~/Documents/Sky Islands/Data/Specimens11092012.csv")
source('~/Documents/Sky Islands/Data/fill_id.R', chdir = TRUE)

spec$essig <- rep("EMEC", nrow(spec))
spec$project <- rep("SI", nrow(spec))
spec$date <- paste("0", spec$date, sep="")

spec$uniqueID <- paste(spec$essig, spec$project, spec$site, spec$date, spec$number, sep='_')

spec$findID <- paste(spec$site, spec$date, spec$number, sep='_')

source('~/Documents/Sky Islands/bomby_morpho.R', chdir = TRUE)

for(i in 1:length(bomby.list)){
	spec <- fill.id(dats=spec, find.id=bomby.list[[1]], sp.id=names(bomby.list)[1], genus.id=NA, family.id="Bombyliidae", determine.date= NA, determiner= NA)
}
