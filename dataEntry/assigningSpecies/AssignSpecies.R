rm(list=ls())
setwd('~/Dropbox/SkyIslands/dataEntry/assigningSpecies')
source('src/AddData.R')

spec <- read.csv(file="../../data/specimenData/cleaned/specimens.csv", as.is=TRUE)
spec$order <- spec$family <- spec$genus <- spec$subgenus <- spec$species <-
  spec$subspecies <-  spec$sex <-  spec$determiner <-
  spec$dateDetermined <-  spec$author <- NA

write.csv(spec,
      file="../../data/specimenData/cleaned/specimens.csv", row.names=FALSE) 

source('Diptera/Syrphidae.R')
add.to.data(sp.ids=sp.ids,
            case='syr', "Syrphidae", no.fam="TRUE", "M. Hauser", "2013")

source('Lepidoptera/Butterflies.R')
add.to.data(sp.ids=sp.ids,
            case='lep', "", no.fam="FALSE", "L.Ponisio", "2013")

source('Hymenoptera/Bees/Andrenidae.R')
add.to.data(sp.ids=sp.ids,
            case='bee', "Andrenidae",  no.fam="TRUE", "J.S. Ascher", "2013")

source('Hymenoptera/Bees/Apidae.R')
add.to.data(sp.ids=sp.ids,
            case='bee', "Apidae", no.fam="TRUE", "J.S. Ascher",  "2013")

source('Hymenoptera/Bees/Bombus.R')
add.to.data(sp.ids=sp.ids,
            case='bee', "Apidae", no.fam="TRUE", "R.W. Thorp", "2013")

source('Hymenoptera/Bees/Halictidae.R')
add.to.data(sp.ids=sp.ids,
            case='bee', "Halictidae", no.fam="TRUE", "J.S. Ascher","2013")

source('Hymenoptera/Bees/Colletidae.R')
add.to.data(sp.ids=sp.ids,
            case='bee', "Colletidae", no.fam="TRUE", "J.S. Ascher","2013")

source('Hymenoptera/Bees/Megachilidae.R')
add.to.data(sp.ids=sp.ids,
            case='bee',"Megachilidae", no.fam="TRUE", "J.S. Ascher","2013")

source('Hummingbirds/Hummingbirds.R')
add.to.data(sp.ids=sp.ids,
            case='humm',"Trochilidae", no.fam="TRUE", "A.J. Rominger","2013")

source('Hymenoptera/Vespids.R')
add.to.data(sp.ids=sp.ids,
            case='was',"Vespidae", no.fam="TRUE", "J. Carpenter","2013")

source('~/Dropbox/SkyIslands/data/relational/src/relational_prep.R')
source('~/Dropbox/SkyIslands/data/relational/src/relational_make.R')
source('~/Dropbox/SkyIslands/data/relational/src/make_traditional.R')
source('~/Dropbox/SkyIslands/analysis/data/dataPrep.R')
