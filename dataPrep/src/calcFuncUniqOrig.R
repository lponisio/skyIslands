library(FD)


## calcFuncUniqOrig <- function(traits, traits.2.keep,
##                              weights,
##                              type = "network",
##                              ...){
##   ## this function calculates trait uniqueness and originality based
##   ## coux et al. FOR A SINGLE COMMUNITY

##   these.traits <- traits[, traits.2.keep]

##   ## convert traits with only 2 catagories to binary, numeric
##   for(this.trait in  traits.2.keep){
##     out.traits <- unique(these.traits[, this.trait])
##     out.traits <- out.traits[!is.na(out.traits)]
##     if(!is.numeric(out.traits)){
##       if(length(out.traits) == 2){
##         print(paste("cat == 2", this.trait))
##         these.traits[, this.trait] <-
##           as.numeric(as.factor(these.traits[,this.trait])) -1
##       } else {
##         print(paste("cat > 2", this.trait))
##         these.traits[, this.trait] <-
##           as.factor(these.traits[,this.trait])
##       }
##     } else {
##       print(paste("numeric", this.trait))
##       these.traits[, this.trait] <-
##         standardize(these.traits[,this.trait])
##     }
##   }
##   ## *********************
##   ## eventually have meadow-level data here, vs. across all SI, can
##   ## save site-level data here as well

##   site.func.mets <- dbFD(these.traits,
##                          w=weights,
##                          corr="lingoes", print.pco=T,
##                           ...)
##   coords <- site.func.mets$x.axes
##   centr <- apply(coords, 2, mean)

##   ## add centroid coords as last row in dataframe
##   coords2 <- data.frame(rbind(coords, centr))
##   rownames(coords2)[dim(coords)[1]+1] <- "centr"

##   ## create a matrix of distances between all species coordinates and
##   ## centroid
##   dists_centr <- as.matrix(dist(coords2, diag=TRUE, upper=TRUE))
##   for (i in 1:dim(dists_centr)[1]) {
##     dists_centr[i,i] <- NA
##   }
##   ## Originality: Distance to centroid of the species present
##   originality <- dists_centr[, dim(dists_centr)[1]]
##   originality <- originality[-(length(originality))]

##   ## Uniqueness: nearest neighbour among present species
##   ## the the minimum disance for each row
##   uniq <- apply(dists_centr[-(dim(dists_centr)[1]),
##                             -(dim(dists_centr)[1])],
##                 1, min, na.rm=TRUE)
##   names(uniq) <- names(originality)
##   out <- cbind(scale(uniq),
##                scale(originality))
##   colnames(out) <- c("uniq", "originality")
##   out <- as.data.frame(out)
##   if(is.null(rownames(these.traits)) | type == "all spec"){
##     out$GenusSpecies <- traits$GenusSpecies
##   } else{
##     out$GenusSpecies <- rownames(traits)
##   }
##   return(out)
## }

calcFuncUniqOrig <- function(traits, traits.2.keep,
                             weights,
                             type = "network",
                             a= NULL,
                             w.abun=TRUE,
                             ...){
  ## this function calculates trait uniqueness and originality based
  ## coux et al. The centroid is calculate for each site, year round. 

  ## drop unwanted traits
  these.traits <- traits[, traits.2.keep]

  ## alphabatize the community matrix
  a <- a[,order(colnames(a))]
  
  ## convert traits with only 2 catagories to binary, numeric
  for(this.trait in  traits.2.keep){
    out.traits <- unique(these.traits[, this.trait])
    out.traits <- out.traits[!is.na(out.traits)]
    if(!is.numeric(out.traits)){
      if(length(out.traits) == 2){
        print(paste("cat == 2", this.trait))
        these.traits[, this.trait] <-
          as.numeric(as.factor(these.traits[,this.trait])) -1
      } else {
        print(paste("cat > 2", this.trait))
        these.traits[, this.trait] <-
          as.factor(these.traits[,this.trait])
      }
    } else {
      print(paste("numeric", this.trait))
      these.traits[, this.trait] <-
        standardize(these.traits[,this.trait])
    }
  }
  ## *********************

  site.func.mets <- dbFD(these.traits,
                         w=weights,
                         corr="lingoes", print.pco=T,
                         a=a, w.abun=w.abun, ...)

  coords <- site.func.mets$x.axes

  coords.comm <- vector(mode="list", length=nrow(a))
  for(i in 1:nrow(a)){
    coords.comm[[i]] <- empty(unlist(c(a[i, ])) * as.matrix(coords))
    ## check
    all(names(unlist(c(a[i, ]))[unlist(c(a[i, ])) > 0]) %in%
        rownames(coords.comm[[i]]))
    all(rownames(coords.comm[[i]]) %in%
        names(unlist(c(a[i, ]))[unlist(c(a[i, ])) > 0]))
    
    
  }

  calcSiteLevelSpMets <- function(coords){
    centr <- apply(coords, 2, mean)
    ## add centroid coords as last row in dataframe
    coords2 <- data.frame(rbind(coords, centr))
    rownames(coords2)[dim(coords)[1]+1] <- "centr"

    ## create a matrix of distances between all species coordinates and
    ## centroid
    dists_centr <- as.matrix(dist(coords2, diag=TRUE, upper=TRUE))
    for (i in 1:dim(dists_centr)[1]) {
      dists_centr[i,i] <- NA
    }
    ## Originality: Distance to centroid of the species present
    originality <- dists_centr[, dim(dists_centr)[1]]
    originality <- originality[-(length(originality))]

    ## Uniqueness: nearest neighbour among present species
    ## the the minimum disance for each row
    uniq <- apply(dists_centr[-(dim(dists_centr)[1]),
                              -(dim(dists_centr)[1])],
                  1, min, na.rm=TRUE)
    names(uniq) <- names(originality)
    out <- cbind(scale(uniq),
                 scale(originality))
    colnames(out) <- c("uniq", "originality")
    out <- as.data.frame(out)
    return(out)
  }

  by.comm.mets <- lapply(coords.comm, calcSiteLevelSpMets)
  names(by.comm.mets) <- rownames(a)

  for(i in 1:length(by.comm.mets)){
    this.name <-  names(by.comm.mets)[i]
    this.comm <- by.comm.mets[[i]]
    this.comm$SiteYearSr <- this.name
    this.comm$GenusSpecies <- rownames(this.comm)
    rownames(this.comm) <- NULL
    by.comm.mets[[i]] <- this.comm
  }

  by.comm.mets <- do.call(rbind, by.comm.mets)
  rownames(by.comm.mets) <- NULL
  
  return(list(fd=   site.func.mets,
              by.comm.mets = by.comm.mets))
}

