
sp.ids <- list(

  Anthidium_maculosum= list(genus="Anthidium",
    subgenus ="", species="maculosum",
    subspecies="", sex="", author = "Cresson",
    temp.id= c("MM_080212_166", "MM_080212_198", "MM_080112_16")),
  
  Megachile_frigida = list(genus="Megachile", subgenus ="",
    species="frigida", subspecies="", sex="m", author = "Smith",
    temp.id=  c("MM_073112_2","JC_071512_156", "JC_071612_89",
      "JC_071612_26")),

   Megachile_relativa = list(genus="Megachile",
    subgenus ="", species="relativa",
    subspecies="", sex="f", author = "Cresson",
    temp.id= c("JC_071712_60", "JC_071712_160")),

  Megachile_policaris = list(genus="Megachile",
    subgenus ="", species="policaris",
    subspecies="", sex="m", author = "Say",
    temp.id= c("PL_081312_95")),

   Megachile_inimica_sayi  = list(genus="Megachile",
    subgenus ="", species="inimica",
    subspecies="sayi", sex="m", author = "Cresson",
    temp.id= c("PL_081112_210")),

   Megachile_melanophaea = list(genus="Megachile",
    subgenus ="", species="melanophaea",
    subspecies="", sex="f", author = "Smith",
    temp.id= c("JC_071712_151", "JC_071612_92", "JC_071712_16",
    "JC_071512_157", "JC_071712_149", "JC_071512_125")),

   Megachile_fidelis = list(genus="Megachile", subgenus ="",
    species="fidelis", subspecies="", sex="f", author = "Cresson",
    temp.id= c("SC_072012_81", "SC_072112_151", "PL_081012_247")),

  Coelioxys_sodalis = list(genus="Coelioxys",
    subgenus ="", species="sodalis",
    subspecies="", sex="f", author = "Cresson",
    temp.id= c("JC_071612_99"))


  )


## ############## read in data and match temp.ids ###################

## lengths <- sapply(sp.ids, function(x) length(x$temp.id))
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
## D$family[ind] <- 'Megachilidae'
## D$genus[ind] <- as.character(dd$Genus)
## D$subgenus[ind] <- as.character(dd$SubGenus)
## D$species[ind] <- as.character(dd$Species)
## D$subspecies[ind] <- as.character(dd$SubSpecies)
## D$sex[ind] <- as.character(dd$Sex)
## D$determiner[ind] <- 'J.S. Ascher'
## D$dateDetermined[ind] <- '2013'

## write.csv(D, file=newfile, row.names=FALSE)
