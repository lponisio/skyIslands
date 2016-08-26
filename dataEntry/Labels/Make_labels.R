
make.labels <- function(morpho.list, spec, table.name, genID){
	## modify date:
	spec$date <- paste("0", spec$date, sep="")

	## create unique ID:
	spec$uniqueID <- paste(rep("EMEC", nrow(spec)),
                       rep("SI", nrow(spec)),
                       spec$site,
                       spec$date,
                       spec$number,
                       sep='_')

	## add columns with the ID to be found in the sourced file
	spec$findID <- paste(spec$site, spec$date, spec$number, sep='_')

	##delist the list of morpho species
	ID.only <- unlist(morpho.list)
	print(length(ID.only))

	##take a subset of the data that is only the morpho species
	sub.dat <- spec[spec$findID %in% ID.only,]

	## merge the site data and the speciment data
	sub.dat <- merge(sub.dat, geo, by="site")
	
	sub.dat <- sub.dat[match(ID.only,sub.dat$findID),]
	
	labels <- sub.dat[, c("genID", "uniqueID", "finalPlantID", "collector", "county", "country", "lat", "long", "state")]
	print(nrow(labels))
	
	write.table(labels, table.name, sep=",", row.names=FALSE)
	
	to.check <- sub.dat[sub.dat$genID != genID, c("findID","genID")]
	
	if(length(ID.only) == nrow(labels)){
		print("TRUE")
		}else{
		print("FALSE")
	}
		
	return(to.check)
}
