add.to.data <- function(sp.ids, case, family, date, data.file) {

    spec.dat <- read.csv(file=data.file, as.is=TRUE)
    lengths <- sapply(sp.ids, function(x) length(x$temp.id))
    cats <- names(sp.ids[[1]])[names(sp.ids[[1]]) != "temp.id"]

    ## make into data-frame
    dd <- sapply(cats, function(cat)
        rep(sapply(sp.ids, function(x) x[[cat]]), lengths))
    TempID <- unlist(sapply(sp.ids, function(x) x$temp.id))

   ## check that no IDs are already present in main data-set
    if(any(TempID %in% spec.dat$SpecimenID[!is.na(spec.dat$Species)])) {
        cat('!!!TempID already present!!!\n')
        dup.temp <- TempID[TempID %in%
                     spec.dat$SpecimenID[!is.na(spec.dat$Species)]]
        print(dup.temp)
        write.csv(bad.ids, file="cleaning/dup_IDs.csv")
    }
    ind <- match(TempID, spec.dat$SpecimenID)

    if(any(is.na(ind))){
      bad.ids <- TempID[is.na(ind)]
      print("bad temp ID")
      print(bad.ids)
        write.csv(bad.ids, file="cleaning/ids_not_in_specimens.csv")
    }

    if(case=='bee') {
        spec.dat$Order[ind] <- 'Hymenoptera'
        spec.dat$GenID[ind] <- 'Bee'
    }
    if(case=='wasp') {
        spec.dat$Order[ind] <- 'Hymenoptera'
        spec.dat$GenID[ind] <- 'Wasp'
    }

    if(case=='lep') {
        spec.dat$Order[ind] <- 'Lepidoptera'
        spec.dat$GenID[ind] <- 'Lep'
    }
    if(case=='beetle') {
        spec.dat$Order[ind] <- 'Coleoptera'
        spec.dat$GenID[ind] <- 'Beetle'
    }
    if(case=='fly') {
        spec.dat$Order[ind] <- 'Diptera'
        spec.dat$GenID[ind] <- 'Fly'
    }

    good.ids <- !is.na(ind)
    ind <- ind[good.ids]
    spec.dat[ind, cats] <- ""
    spec.dat[ind, cats] <- dd[good.ids, cats]
    spec.dat$DateDetermined[ind] <- date

    if(case != "beetle"){
        spec.dat$Family[ind] <- family
    }
    write.csv(spec.dat,
              file=data.file,
              row.names=FALSE)
}



