
sp.ids <- list(

  Hylaeus_annulatus = list(genus="Hylaeus", subgenus ="",
    species="annulatus", subspecies="", sex="f", author = "Linnaeus",
    temp.id= c("JC_071712_158", "JC_071612_23", "JC_071212_52",
    "JC_071712_92", "JC_071612_125", "JC_071712_177", "SC_072012_289",
    "SC_072712_171", "JC_071212_47","SC_072812_207",
    "SC_072812_210")),

   Colletes_skinneri= list(genus="Colletes", subgenus ="",
    species="skinneri", subspecies="", sex="f", author = "Viereck",
    temp.id= c("PL_081312_133","MM_080412_4","MM_080212_240")),

    Colletes_gilensis= list(genus="Colletes", subgenus ="",
    species="gilensis", subspecies="", sex="f", author = "Cockerell",
    temp.id= c("SC_072712_58"))

  )

############## read in data and match temp.ids ###################

##lengths <- sapply(sp.ids, function(x) length(x$temp.id))
## Genus <- sapply(sp.ids, function(x) x$genus)
## SubGenus <- sapply(sp.ids, function(x) x$subgenus)
## Species <- sapply(sp.ids, function(x) x$species)
## SubSpecies <- sapply(sp.ids, function(x) x$subspecies)
## Sex <- sapply(sp.ids, function(x) x$sex)
## TempID <- sapply(sp.ids, function(x) x$temp.id)

## dd <- data.frame(TempID=unlist(TempID),
##                  Genus=rep(Genus, lengths),
##                  SubGenus=rep(SubGenus, lengths),
##                  Species=rep(Species, lengths),
##                  SubSpecies=rep(SubSpecies, lengths),
##                  Sex=rep(Sex, lengths))

## ## lood data:
## D <- read.csv(file=newfile, as.is=TRUE)


## ind <- match(dd$TempID, D$temp.id)
## D$order[ind] <- 'Hymenoptera'
## D$family[ind] <- 'Colletidae'
## D$genus[ind] <- as.character(dd$Genus)
## D$subgenus[ind] <- as.character(dd$SubGenus)
## D$species[ind] <- as.character(dd$Species)
## D$subspecies[ind] <- as.character(dd$SubSpecies)
## D$sex[ind] <- as.character(dd$Sex)
## D$determiner[ind] <- 'J.S. Ascher'
## D$dateDetermined[ind] <- '2013'

## write.csv(D, file=newfile, row.names=FALSE)
