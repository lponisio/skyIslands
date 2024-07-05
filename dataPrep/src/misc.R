
## Function for updating plant names based on output of Taxonstand

## fixPlantNames <- function(data.to.fix, ## specimen or veg data
##                           col.name.to.fix="PlantGenusSpecies", ## column name of plant species
##                           checked.plant.names ## output of Taxonstand
##                           ){
##     plant.names <- data.frame(newNames = fix.white.space(paste(checked.plant.names$New.Genus,
##                                                                checked.plant.names$New.Species,
##                                                                checked.plant.names$New.Infraspecific.rank,
##                                                                checked.plant.names$New.Infraspecific)),
##                               oldNames = checked.plant.names$Taxon)
##     data.to.fix$PlantFamily  <- NA
##     data.to.fix$PlantFamily <- checked.plant.names$Family[match(data.to.fix[, col.name.to.fix],
##                                                                 checked.plant.names$Taxon)]

##     data.to.fix$PlantGenus <- NA
##     data.to.fix$PlantGenus <- checked.plant.names$New.Genus[match(data.to.fix[, col.name.to.fix],
##                                                                   checked.plant.names$Taxon)]
    
##     data.to.fix[, col.name.to.fix] <-
##         plant.names$newNames[match(data.to.fix[, col.name.to.fix],
##                                    plant.names$oldNames)]
##     return(data.to.fix)

## }



fixPlantNamesgBIF <- function(data.to.fix, ## specimen or veg data
                          col.name.to.fix="PlantGenusSpecies", ## column name of plant species
                          checked.plant.names ## output of Taxonstand
                          ){
    plant.names <- data.frame(newNames = checked.plant.names$canonicalName,
                              oldNames = checked.plant.names$verbatimScientificName)
    data.to.fix$PlantFamily  <- NA
    data.to.fix$PlantFamily <- checked.plant.names$family[match(data.to.fix[, col.name.to.fix],
                                                   checked.plant.names$verbatimScientificName)]

    data.to.fix$PlantGenus <- NA
    data.to.fix$PlantGenus <- checked.plant.names$genus[match(data.to.fix[, col.name.to.fix],
                                                   checked.plant.names$verbatimScientificName)]

    data.to.fix[, col.name.to.fix] <-
        checked.plant.names$canonicalName[
            match(data.to.fix[, col.name.to.fix],
                  checked.plant.names$verbatimScientificName)]
   
    return(data.to.fix)

}


## standardize a vector
standardize <- function(x)
  (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)


## This functions takes site-species-abundance data and creates a
## matrix where the sites are columns and the rows are species.

samp2site.spp <- function(site, spp, abund, FUN=mean) {
  x <- tapply(abund, list(site = site, spp = spp), FUN)
  x[is.na(x)] <- 0
  return(x)
}

## does the reverse of samp2site
comm.mat2sample <-  function (z) {
  temp <- data.frame(expand.grid(dimnames(z))[1:2],
                     as.vector(as.matrix(z)))
  temp <- temp[sort.list(temp[, 1]), ]
  data.frame(Site = temp[, 1], Samp = temp[, 3],
             Date = temp[, 2])
}


## function to clean up white-space in a column of data (replaces all
## instances of white-space with " " and empty cells with ""
fix.white.space <- function(d) {
  d <- as.character(d)
  remove.first <- function(s) substr(s, 2, nchar(s))
  d <- gsub("      ", " ", d, fixed=TRUE)
  d <- gsub("     ", " ", d, fixed=TRUE)
  d <- gsub("    ", " ", d, fixed=TRUE)
  d <- gsub("   ", " ", d, fixed=TRUE)
  d <- gsub("  ", " ", d, fixed=TRUE)

  tmp <- strsplit(as.character(d), " ")
  d <- sapply(tmp, function(x) paste(x, collapse=" "))

  first <- substr(d, 1, 1)
  d[first==" "] <- remove.first(d[first==" "])
  d
}


## return sorted unique values
id <- function(x) unique(sort(x))

## function to make pollinator visitation matrices
make.mats <- function(pollinator.id,
                      null.mat,
                      pollinator,
                      var1,
                      var2,
                      occ) {
  make.mat <- function(P) {
    var1 <- var1[pollinator==P]
    var2 <- var2[pollinator==P]
    occ <- occ[pollinator==P]
    null.mat[unique(var1),
             unique(as.character(var2))][!is.na(null.mat)] <- occ
    null.mat
  }
  mats <- lapply(pollinator.id, function(x) make.mat(x))
  names(mats) <- pollinator.id
  mats
}


catchDups <- function(indiv.comm){
    ## this function combine columns with duplicate names in the
    ## community matrix
    dups  <-
        na.omit(unique(colnames(indiv.comm)[(duplicated(colnames(indiv.comm)))]))
    #browser()
    print(length(dups))
    for(dup in dups){
        columnstomerge <- indiv.comm[,colnames(indiv.comm) == dup]
        
        indiv.comm <- indiv.comm[,!colnames(indiv.comm) == dup]
        indiv.comm <- cbind(indiv.comm, apply(columnstomerge, 1, sum))
        colnames(indiv.comm)[colnames(indiv.comm) == ""] <- dup
    }
    return(indiv.comm)
}


makeComm <- function(taxonomy, features, feature.col="Feature.ID"){
    ## this function make a specimne by species matrix
    rownames(taxonomy) <- features$Taxon[match(rownames(taxonomy),
                                               features[, feature.col])]
    
    
    indiv.comm <- t(taxonomy)
    #browser()
    return(indiv.comm)
}
