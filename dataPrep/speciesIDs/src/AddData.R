add.to.data <- function(sp.ids, case, family, date, data.file) {
    D <- read.csv(file=data.file, as.is=TRUE)

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

    ## check that no IDs are already present in main data-set
    if(any(TempID %in% D$UniqueID[!is.na(D$Species)])) {
        cat('!!!TempID already present!!!\n')
        print(TempID[TempID %in% D$Unique[!is.na(D$Species)]])
    }
    ind <- match(TempID, D$UniqueID)

    if(any(is.na(ind))){
        print(paste("bad temp ID", TempID[is.na(ind)]))
    }
    if(case=='bee') {
        D$Order[ind] <- 'Hymenoptera'
        D$GenID[ind] <- 'Bee'
    }
    if(case=='wasp') {
        D$Order[ind] <- 'Hymenoptera'
        D$GenID[ind] <- 'Wasp'
    }

    if(case=='lep') {
        D$Order[ind] <- 'Lepidoptera'
        D$GenID[ind] <- 'Lep'
    }
    if(case=='beetle') {
        D$Order[ind] <- 'Coleoptera'
        D$GenID[ind] <- 'Beetle'
    }
    if(case=='fly') {
        D$Order[ind] <- 'Diptera'
        D$GenID[ind] <- 'Fly'
    }
    D[ind,cats] <- dd[,cats]
    D$DateDetermined[ind] <- date

    if(case != "beetle"){
        D$Family[ind] <- family
    }
    write.csv(D,
              file=data.file,
              row.names=FALSE)
}



