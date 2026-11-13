## Wasps occur in several families even when their IDs were archived in
## Vespids.R or Sphecidae.R. Do not infer Family from the file name alone.
## This is a deliberately explicit, auditable lookup for the wasp genera
## represented in the 2012, 2017 and 2018 ID scripts.
## Crabronidae is used sensu lato (consistent with the existing ID tables).
wasp_family_by_genus <- c(
  Ancistrocerus='Vespidae', Dolichovespula='Vespidae',
  Eumenes='Vespidae', Euodynerus='Vespidae',
  Odynerus='Vespidae', Polistes='Vespidae',
  Pseudepipona='Vespidae', Stenodynerus='Vespidae',
  Pepsis='Pompilidae', Cryptocheilus='Pompilidae',
  Chirodamus='Pompilidae',
  Scolia='Scoliidae', Triscolia='Scoliidae',
  Ammophila='Sphecidae', Podalonia='Sphecidae', Sphex='Sphecidae',
  Bembix='Crabronidae', Crabro='Crabronidae',
  Ectemnius='Crabronidae', Philanthus='Crabronidae',
  Tachytes='Crabronidae',
  ## Family-only IDs in 2018/Vespids.R are stored in the Genus column.
  Crabronidae='Crabronidae', Scoliidae='Scoliidae'
)

add.to.data <- function(sp.ids, case, family, date, data.file) {
  spec.dat <- read.csv(file=data.file, as.is=TRUE)
  lengths <- sapply(sp.ids, function(x) length(x$temp.id))
  cats <- names(sp.ids[[1]])[names(sp.ids[[1]]) != "temp.id"]

  ## make into data-frame
  dd <- sapply(cats, function(cat)
    rep(sapply(sp.ids, function(x) x[[cat]]), lengths))
  TempID <- unlist(sapply(sp.ids, function(x) x$temp.id))
  if(length(TempID == 1)){
    names(dd) <- cats
  }

  ## check for duplicates
  if(sum(table(TempID) > 1) > 0) {
    cat('!!!Duplicate TempID found!!!\n')
    print(family)
    tab.dup <- table(TempID)[table(TempID) > 1]
    print(tab.dup)
    write.csv(tab.dup, file=sprintf(
                         "cleaning/%s_%s_%s_dup_IDs_same.csv",
                         case,
                         family,
                         date))
  }

  ## check that no IDs are already present in main data-set
  if(any(TempID %in% spec.dat$UniqueID[!is.na(spec.dat$Species)])) {
    cat('!!!TempID already present!!!\n')
    dup.temp <- TempID[TempID %in%
                       spec.dat$UniqueID[!is.na(spec.dat$Species)]]
    print(dup.temp)
    write.csv(dup.temp, file=sprintf(
                          "cleaning/%s_%s_%s_dup_IDs_diff.csv",
                          case,
                          family,
                          date))
  }
  ind <- match(TempID, spec.dat$UniqueID)

  if(any(is.na(ind))){
    bad.ids <- TempID[is.na(ind)]
    print("bad temp ID")
    print(bad.ids)
    write.csv(bad.ids, file=sprintf(
                         "cleaning/%s_%s_%s_ids_not_in_raw.csv",
                         case,
                         family,
                         date))
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
  if(length(good.ids) > 1){
    spec.dat[ind, cats] <- dd[good.ids, cats]
  } else{
    spec.dat[ind, cats] <- dd[cats]
  }
  spec.dat$DateDetermined[ind] <- date

  ## For Lepidoptera, the ID lists give Family per taxon but callers pass
  ## family="". Do not overwrite those species-level family IDs with blanks:
  ## the downstream butterfly / hawk-moth guilds are family-defined.
  if(case != "beetle" && !is.na(family) && nzchar(family)){
    spec.dat$Family[ind] <- family
  }

  ## Family assignment for wasps must use the specimen's genus, not the
  ## broad file label ('Vespidae' or 'Sphecidae') supplied by AssignSpecies.
  ## Fail loudly on new / unmapped genera rather than writing a false Family.
  if (case == 'wasp' && length(ind)) {
    wasp.genera <- trimws(as.character(spec.dat$Genus[ind]))
    mapped.families <- unname(wasp_family_by_genus[wasp.genera])
    if (anyNA(mapped.families)) {
      stop('Unmapped wasp genus in AddData.R: ',
           paste(unique(wasp.genera[is.na(mapped.families)]),
                 collapse=', '),
           '. Add the verified family to wasp_family_by_genus.')
    }
    spec.dat$Family[ind] <- mapped.families
  }
  write.csv(spec.dat,
            file=data.file,
            row.names=FALSE)
}



