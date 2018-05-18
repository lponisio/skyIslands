
sp.ids <- list(

  Anthidium_maculosum= list(Genus="Anthidium",
    SubGenus ="", Species="maculosum",
    SubSpecies="", Sex="", Author = "Cresson",
    temp.id= c("MM_080212_166", "MM_080212_198", "MM_080112_16")),
  
  Megachile_frigida = list(Genus="Megachile", SubGenus ="",
    Species="frigida", SubSpecies="", Sex="m", Author = "Smith",
    temp.id=  c("MM_073112_2","JC_071512_156", "JC_071612_89",
      "JC_071612_26")),

   Megachile_relativa = list(Genus="Megachile",
    SubGenus ="", Species="relativa",
    SubSpecies="", Sex="f", Author = "Cresson",
    temp.id= c("JC_071712_60", "JC_071712_160")),

  Megachile_policaris = list(Genus="Megachile",
    SubGenus ="", Species="policaris",
    SubSpecies="", Sex="m", Author = "Say",
    temp.id= c("PL_081312_95")),

   Megachile_inimica_sayi  = list(Genus="Megachile",
    SubGenus ="", Species="inimica",
    SubSpecies="sayi", Sex="m", Author = "Cresson",
    temp.id= c("PL_081112_210")),

   Megachile_melanophaea = list(Genus="Megachile",
    SubGenus ="", Species="melanophaea",
    SubSpecies="", Sex="f", Author = "Smith",
    temp.id= c("JC_071712_151", "JC_071612_92", "JC_071712_16",
    "JC_071512_157", "JC_071712_149", "JC_071512_125")),

   Megachile_fidelis = list(Genus="Megachile", SubGenus ="",
    Species="fidelis", SubSpecies="", Sex="f", Author = "Cresson",
    temp.id= c("SC_072012_81", "SC_072112_151", "PL_081012_247")),

  Coelioxys_sodalis = list(Genus="Coelioxys",
    SubGenus ="", Species="sodalis",
    SubSpecies="", Sex="f", Author = "Cresson",
    temp.id= c("JC_071612_99"))


  )


## ############## read in data and match temp.ids ###################

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
## D$order[ind] <- 'Hymenoptera'
## D$family[ind] <- 'Megachilidae'
## D$Genus[ind] <- as.character(dd$Genus)
## D$SubGenus[ind] <- as.character(dd$SubGenus)
## D$Species[ind] <- as.character(dd$Species)
## D$SubSpecies[ind] <- as.character(dd$SubSpecies)
## D$Sex[ind] <- as.character(dd$Sex)
## D$Determiner[ind] <- 'J.S. Ascher'
## D$dateDetermined[ind] <- '2013'

## write.csv(D, file=newfile, row.names=FALSE)
