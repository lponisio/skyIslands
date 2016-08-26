add.to.data <- function(sp.ids, case, family, no.fam, determiner, date) {
  lengths <- sapply(sp.ids, function(x) length(x$temp.id))
  cats <- names(sp.ids[[1]])[-length(sp.ids[[1]])]

  ## make into data-frame
  dd <- sapply(cats, function(cat)
               rep(sapply(sp.ids, function(x) x[[cat]]), lengths))
  TempID <- unlist(sapply(sp.ids, function(x) x$temp.id))
  
  ## check for duplicates
  if(sum(table(TempID)>1)>0) {
    cat('!!!Duplicate TempID found!!!\n')
    print(table(TempID)[table(TempID)>1])
  }

  ## lood data:
  D <- read.csv(file="~/Dropbox/SkyIslands/data/specimenData/cleaned/specimens.csv", as.is=TRUE)
  ## check that no IDs are already present in main data-set
  if(any(TempID %in% D$TempID[!is.na(D$species)])) {
    cat('!!!TempID already present!!!\n')
    print(TempID[TempID %in% D$TempID[!is.na(D$species)]])
  }
  ind <- match(TempID, D$temp.id)
  if(case=='bee') {
    D$order[ind] <- 'Hymenoptera'
  }
  if(case=='was') {
    D$order[ind] <- 'Hymenoptera'
  }
  if(case=='syr') {
    D$order[ind] <- 'Diptera'
    D$family[ind] <- 'Syrphidae'
  }
  if(case=='lep') {
    D$order[ind] <- 'Lepidoptera'
  }
  if(case=="humm"){
    D$order[ind] <- 'Apodiformes'
    D$family[ind] <- 'Trochilidae'
  }
  D[ind,cats] <- dd[,cats]
  if(no.fam){
    D$family[ind] <- family
  }
  D$determiner[ind] <- determiner
  D$dateDetermined[ind] <- date

  write.csv(D,
            file="~/Dropbox/SkyIslands/data/specimenData/cleaned/specimens.csv", row.names=FALSE) }

