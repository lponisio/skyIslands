setwd('~/Dropbox/skyislands/dataEntry/speciesIDs')
source('src/AddData.R')

spec.data.file <- "../../data/relational/data/original/specimens.csv"

spec <- read.csv(file=spec.data.file, as.is=TRUE)
spec$order <- spec$family <- spec$genus <- spec$subgenus <- spec$species <-
  spec$subspecies <-  spec$sex <-  spec$determiner <-
  spec$dateDetermined <-  spec$author <- NA


## 2017 field season
source('IDs/2017/Apidae.R')
add.to.data(sp.ids=sp.ids,
            case='bee', "Apidae","2018",
            D=spec,
            data.file=spec.data.file)

source('IDs/2017/Colletidae.R')
add.to.data(sp.ids=sp.ids,
            case='bee', "Colletidae","2018",
            D=spec,
            data.file=spec.data.file)

source('IDs/2017/Halictidae.R')
add.to.data(sp.ids=sp.ids,
            case='bee', "Halictidae","2018",
            D=spec,
            data.file=spec.data.file)

source('IDs/2017/Megachilidae.R')
add.to.data(sp.ids=sp.ids,
            case='bee', "Megachilidae","2018", D=spec,
            data.file=spec.data.file)

## 2012 field season
source('IDs/2012/Andrenidae.R')
add.to.data(sp.ids=sp.ids,
            case='bee', "Andrenidae","2013-2015",
             D=spec,
            data.file=spec.data.file)

source('IDs/2012/Apidae.R')
add.to.data(sp.ids=sp.ids,
            case='bee', "Apidae","2013-2015",
             D=spec,
            data.file=spec.data.file)

source('IDs/2012/Colletidae.R')
add.to.data(sp.ids=sp.ids,
            case='bee', "Colletidae","2013-2015",
             D=spec,
            data.file=spec.data.file)

source('IDs/2012/Halictidae.R')
add.to.data(sp.ids=sp.ids,
            case='bee', "Halictidae","2013-2015",
             D=spec,
            data.file=spec.data.file)

source('IDs/2012/Megachilidae.R')
add.to.data(sp.ids=sp.ids,
            case='bee',"Megachilidae","2013-2015",
             D=spec,
            data.file=spec.data.file)

source('IDs/2012/Vespids.R')
add.to.data(sp.ids=sp.ids,
            case='wasp',"Vespidae","2015",
             D=spec,
            data.file=spec.data.file)

source('IDs/2012/Butterflies.R')
add.to.data(sp.ids=sp.ids,
            case='lep',"","2012",
             D=spec,
            data.file=spec.data.file)

source('IDs/2012/Syrphidae.R')
add.to.data(sp.ids=sp.ids,
            case='fly',"Syrphidae","2012",
             D=spec,
            data.file=spec.data.file)
