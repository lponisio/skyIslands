
save.dir <- "../../skyIslands/dataPrep/speciesIDs"

source(file.path(save.dir, 'src/AddData.R'))


spec.data.file <- "relational/original/specimens.csv"

spec <- read.csv(file=spec.data.file, as.is=TRUE)

spec$Order <- spec$Family <- spec$Genus <- spec$SubGenus <- spec$Species <-
  spec$SubSpecies <-  spec$Sex <-  spec$Determiner <-
      spec$DateDetermined <-  spec$Author <- NA

write.csv(spec, file=spec.data.file, row.names=FALSE)




## 2022 field season
## source(file.path(save.dir, 'IDs/2021/Syrphidae.R'))
## add.to.data(sp.ids=sp.ids,
##             case='fly',"Syrphidae", "2021",
##             data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2022/Apidae.R'))
add.to.data(sp.ids=sp.ids, case='bee', "Apidae","2022",
data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2022/Halictidae.R'))
add.to.data(sp.ids=sp.ids,
            case='bee', "Halictidae","2022",
            data.file=spec.data.file)

## source(file.path(save.dir, 'IDs/2022/Megachilidae.R'))
## add.to.data(sp.ids=sp.ids,
##             case='bee', "Megachilidae","2022",
##             data.file=spec.data.file)

## source(file.path(save.dir, 'IDs/2022/Colletidae.R'))
## add.to.data(sp.ids=sp.ids,
##             case='bee', "Colletidae","2022",
##             data.file=spec.data.file)

## source(file.path(save.dir, 'IDs/2022/Sphecidae.R'))
## add.to.data(sp.ids=sp.ids,
##             case='wasp', "Sphecidae","2022",
##             data.file=spec.data.file)



## 2021 field season
## source(file.path(save.dir, 'IDs/2021/Syrphidae.R'))
## add.to.data(sp.ids=sp.ids,
##             case='fly',"Syrphidae", "2021",
##             data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2021/Apidae.R'))
add.to.data(sp.ids=sp.ids, case='bee', "Apidae","2021",
data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2021/Halictidae.R'))
add.to.data(sp.ids=sp.ids,
            case='bee', "Halictidae","2021",
            data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2021/Megachilidae.R'))
add.to.data(sp.ids=sp.ids,
            case='bee', "Megachilidae","2021",
            data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2021/Colletidae.R'))
add.to.data(sp.ids=sp.ids,
            case='bee', "Colletidae","2021",
            data.file=spec.data.file)


## source(file.path(save.dir, 'IDs/2021/Sphecidae.R'))
## add.to.data(sp.ids=sp.ids,
##             case='wasp', "Sphecidae","2021",
##             data.file=spec.data.file)




## 2018 field season
source(file.path(save.dir, 'IDs/2018/Syrphidae.R'))
add.to.data(sp.ids=sp.ids,
            case='fly',"Syrphidae", "2019",
            data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2018/Apidae.R'))
add.to.data(sp.ids=sp.ids, case='bee', "Apidae","2019",
data.file=spec.data.file)


source(file.path(save.dir, 'IDs/2018/Colletidae.R'))
add.to.data(sp.ids=sp.ids, case='bee', "Colletidae","2019",
data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2018/Halictidae.R'))
add.to.data(sp.ids=sp.ids,
            case='bee', "Halictidae","2019",
            data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2018/Megachilidae.R'))
add.to.data(sp.ids=sp.ids,
            case='bee', "Megachilidae","2019",
            data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2018/Sphecidae.R'))
add.to.data(sp.ids=sp.ids,
            case='wasp', "Sphecidae","2020",
            data.file=spec.data.file)


source(file.path(save.dir, 'IDs/2018/Bombyliidae.R'))
add.to.data(sp.ids=sp.ids,
            case='fly', "Bombyliidae","2020",
            data.file=spec.data.file)



## 2017 field season
source(file.path(save.dir, 'IDs/2017/Syrphidae.R'))
add.to.data(sp.ids=sp.ids,
            case='fly',"Syrphidae", "2019",
            data.file=spec.data.file)


source(file.path(save.dir, 'IDs/2017/Apidae.R'))
add.to.data(sp.ids=sp.ids,
            case='bee', "Apidae","2018",
            data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2017/Colletidae.R'))
add.to.data(sp.ids=sp.ids,
            case='bee', "Colletidae","2018",
            data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2017/Halictidae.R'))
add.to.data(sp.ids=sp.ids,
            case='bee', "Halictidae","2018",
            data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2017/Megachilidae.R'))
add.to.data(sp.ids=sp.ids,
            case='bee', "Megachilidae","2018",
            data.file=spec.data.file)


source(file.path(save.dir, 'IDs/2017/Sphecidae.R'))
add.to.data(sp.ids=sp.ids,
            case='wasp', "Sphecidae","2020",
            data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2017/Bombyliidae.R'))
add.to.data(sp.ids=sp.ids,
            case='fly', "Bombyliidae","2020",
            data.file=spec.data.file)



## 2012 field season
source(file.path(save.dir, 'IDs/2012/Apidae.R'))
add.to.data(sp.ids=sp.ids,
            case='bee', "Apidae","2013-2015",
            data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2012/Colletidae.R'))
add.to.data(sp.ids=sp.ids,
            case='bee', "Colletidae","2013-2015",
            data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2012/Halictidae.R'))
add.to.data(sp.ids=sp.ids,
            case='bee', "Halictidae","2013-2015",
            data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2012/Megachilidae.R'))
add.to.data(sp.ids=sp.ids,
            case='bee',"Megachilidae","2013-2015",
            data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2012/Vespids.R'))
add.to.data(sp.ids=sp.ids,
            case='wasp',"Vespidae","2015",
            data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2012/Butterflies.R'))
add.to.data(sp.ids=sp.ids,
            case='lep',"","2012",
            data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2012/Syrphidae.R'))
add.to.data(sp.ids=sp.ids,
            case='fly',"Syrphidae", "2012",
            data.file=spec.data.file)


source(file.path(save.dir, 'IDs/2012/Sphecidae.R'))
add.to.data(sp.ids=sp.ids,
            case='wasp', "Sphecidae","2020",
            data.file=spec.data.file)

source(file.path(save.dir, 'IDs/2012/Bombyliidae.R'))
add.to.data(sp.ids=sp.ids,
            case='fly', "Bombyliidae","2020",
            data.file=spec.data.file)


