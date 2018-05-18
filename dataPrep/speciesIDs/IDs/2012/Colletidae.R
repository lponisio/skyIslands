
sp.ids <- list(

  Hylaeus_annulatus = list(Genus="Hylaeus", SubGenus ="",
    Species="annulatus", SubSpecies="", Sex="f", Author = "Linnaeus",
    temp.id= c("JC_071712_158", "JC_071612_23", "JC_071212_52",
    "JC_071712_92", "JC_071612_125", "JC_071712_177", "SC_072012_289",
    "SC_072712_171", "JC_071212_47","SC_072812_207",
    "SC_072812_210")),

   Colletes_skinneri= list(Genus="Colletes", SubGenus ="",
    Species="skinneri", SubSpecies="", Sex="f", Author = "Viereck",
    temp.id= c("PL_081312_133","MM_080412_4","MM_080212_240")),

    Colletes_gilensis= list(Genus="Colletes", SubGenus ="",
    Species="gilensis", SubSpecies="", Sex="f", Author = "Cockerell",
    temp.id= c("SC_072712_58"))

  )

############## read in data and match temp.ids ###################

##lengths <- sapply(sp.ids, function(x) length(x$temp.id))
## Genus <- sapply(sp.ids, function(x) x$Genus)
## SubGenus <- sapply(sp.ids, function(x) x$SubGenus)
## Species <- sapply(sp.ids, function(x) x$Species)
## SubSpecies <- sapply(sp.ids, function(x) x$SubSpecies)
## Sex <- sapply(sp.ids, function(x) x$Sex)
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
## D$Genus[ind] <- as.character(dd$Genus)
## D$SubGenus[ind] <- as.character(dd$SubGenus)
## D$Species[ind] <- as.character(dd$Species)
## D$SubSpecies[ind] <- as.character(dd$SubSpecies)
## D$Sex[ind] <- as.character(dd$Sex)
## D$Determiner[ind] <- 'J.S. Ascher'
## D$dateDetermined[ind] <- '2013'

## write.csv(D, file=newfile, row.names=FALSE)
