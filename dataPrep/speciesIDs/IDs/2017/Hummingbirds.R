
sp.ids <- list(

  Selasphorus_platycercus= list(Genus="Selasphorus", SubGenus ="",
    Species="platycercus", SubSpecies="", Sex="", Author = "",
    temp.id= c("JC_071212_1001", "JC_071412_1002", "JC_071412_1003",
    "JC_071412_1004", "JC_071412_1005", "JC_071412_1006",
    "JC_071412_1007", "JC_071512_1008", "JC_071512_1009",
    "JC_071512_1010", "JC_071512_1011", "JC_071512_1012",
    "JC_071512_1013", "JC_071612_1014", "JC_071612_1015",
    "JC_071612_1016", "JC_071612_1017", "JC_071612_1018",
    "JC_071612_1019", "JC_071712_1020", "JC_071712_1021",
    "JC_071712_1022", "JC_071712_1023", "JC_071712_1024",
    "JC_071712_1025", "SC_072012_1026", "SC_072012_1027",
    "SC_072012_1028", "SC_072012_1029", "SC_072112_1030",
    "SC_072112_1031", "SC_072112_1032", "CH_082112_1033",
    "CH_082112_1034", "CH_082212_1035", "CH_082212_1036")),

  Selasphorus_rufus= list(Genus="Selasphorus", SubGenus ="",
    Species="rufus", SubSpecies="", Sex="", Author = "", temp.id=
    c("CH_082112_1037", "CH_082112_1038", "CH_082212_1039",
    "CH_082212_1040")))

## lengths <- sapply(sp.ids, function(x) length(x$temp.id))
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
## D$order[ind] <- 'Apodiformes'
## D$family[ind] <- 'Trochilidae'
## D$Genus[ind] <- as.character(dd$Genus)
## D$SubGenus[ind] <- as.character(dd$SubGenus)
## D$Species[ind] <- as.character(dd$Species)
## D$SubSpecies[ind] <- as.character(dd$SubSpecies)
## D$Sex[ind] <- as.character(dd$Sex)
## D$determiner[ind] <- 'A.J. Rominger'
## D$dateDetermined[ind] <- '2013'

## write.csv(D, file=newfile, row.names=FALSE)
