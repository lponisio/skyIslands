rm(list=ls())
setwd("~/Dropbox/SI_data_entry")
source("Make_labels.R")
spec <- read.csv("Specimens03122013.csv")
geo <- read.csv("Geography.csv")


####################################

##morpho.list is the name of the list of morphospecies
##table.name is the name of the table to be written, don't forget the .csv
## genID is the general id of the specimens to be labeled, i.d. "wasp", "syrphid" 

#vespids

source("Species_ids/vespidae_morpho.R")

make.labels(morpho.list=morpho.list, spec=spec, table.name="vespid.labels.csv", genID="wasp")

#syrphids 

source("Species_ids/syrphid_morpho.R")

make.labels(morpho.list=syrphid.list, spec=spec, table.name="syrphid.labels.csv", genID="syrphid")

#bees
source('Species_ids/Bee_morpho.R')

make.labels(morpho.list=bee.morpho, spec=spec, table.name="bee.labels.csv", genID="bee")

#bombus
source('Species_ids/Bombus_morpho.R')
make.labels(morpho.list=bombus.morpho, spec=spec, table.name="bombus.labels.csv", genID="bumble bee")

##redo bombus labels

source('Species_ids/Bombus_redo.R')
make.labels(morpho.list=redo.bomb, spec=spec, table.name="bombus.labels.csv", genID="bumble bee")

##honeybees
source('Species_ids/honey_bees.R')
make.labels(morpho.list=honeybees, spec=spec, table.name="honeybee_labels.csv", genID="honey bee")

##syrphids
source('Species_ids/syrphid_morpho.R')
make.labels(morpho.list=syrphid.list, spec=spec, table.name="syrphid_labels.csv", genID="syrphid")

################################################

#for Kelly
source('Species_ids/redo_ves.R')

make.labels(morpho.list=redo.vespid, spec=spec, table.name="redo_vespids.labels.csv", genID="wasp")


## easy bees 1
source('Species_ids/Easy_bees.R')
make.labels(morpho.list=easy.bees, spec=spec, table.name="agopostemen_anthophora_xylacopa.csv", genID="bee")


### megechile and others
source('Species_ids/Megachile_morpho.R')
make.labels(morpho.list=megachile.morpho, spec=spec, table.name="megachile_others.csv", genID="bee")

##melissodes and others

source('Species_ids/melissodes_sphecodes.R')
make.labels(morpho.list=melissodes, spec=spec, table.name="melissodes_others.csv", genID="bee")

## bombyliids

source('Species_ids/bomby_morpho.R')
make.labels(morpho.list=bomby.list, spec=spec, table.name="bomby_labels.csv", genID="bombyliid")


##sphecids and crabronids
source('Species_ids/Crabronidae_morpho.R')
make.labels(morpho.list=wasp.morpho, spec=spec, table.name="wasp_labels.csv", genID="wasp")

## andrena

source('Species_ids/andrena_morpho.R')
make.labels(morpho.list=mostly.andrena, spec=spec, table.name="andrena_labels.csv", genID="small dark bee")

### butterflies
source('Species_ids/Butterfly_morpho (1).R')
make.labels(morpho.list=butterfly.list, spec=spec, table.name="butterfly_labels.csv", genID="butterfly")
butterfly.spec <- read.csv("butterfly_labels.csv")


##
source('Species_ids/small_bees.R')
make.labels(morpho.list=small.bees, spec=spec, table.name="smallBee_labels.csv", genID="small dark bee")


##ichnumomids, brachonids, pompilids, ants
source('Species_ids/ich_brach_beetle.R')
make.labels(morpho.list=misc.list, spec=spec, table.name="IchBrachPomp_labels.csv", genID="wasp")

##misc flies!
source('Species_ids/Tachinids_morpho_4_16.R')

make.labels(morpho.list=flies.list, spec=spec, table.name="fly_labels.csv", genID="fly")

##bee checking
source('Final species ids/Apidae.R')
make.labels(morpho.list=apidae, spec=spec, table.name="apidae.csv", genID="bee")

##bomby checking
source('Species_ids/bomby_morpho_04_09.R')
make.labels(morpho.list=bomby.list, spec=spec, table.name="bomby_labels.csv", genID="bombyliid")


source("Final species ids/Coleoptera/Beetles.R")
make.labels(morpho.list=beetles, spec=spec, table.name="beetles_labels.csv", genID="beetle")

source('Species_ids/missing.R', chdir = TRUE)
make.labels(morpho.list=missing, spec=spec, table.name="missing_labels.csv", genID="bee")

source('Species_ids/Bombus_morpho_050213.R', chdir = TRUE)
make.labels(bombus.morpho, spec=spec, table.name="bombus_labels.csv", genID="bumble bee")

source('Final species ids/Lepidoptera/Moths.R')
make.labels(moths.morpho, spec=spec, table.name="moth_labels.csv", genID="moth")


source('Final species ids/Vespoidea/Scoliidae.R')
make.labels(scolia, spec=spec, table.name="scolia_labels.csv", genID="wasp")

source('Final species ids/Apoidea/Crabronidae.R')

make.labels(wasp.morpho, spec=spec, table.name="crabronid_labels.csv", genID="wasp")


source('Final species ids/Apoidea/Sphecidae.R', chdir = TRUE)

make.labels(sphecidae, spec=spec, table.name="specid_labels.csv", genID="wasp")

source('Final species ids/Hymenoptera/Vespoidea/Vespidae.R')
make.labels(vespidae_morpho, spec=spec, table.name="ves_labels.csv", genID="wasp")


source('Final species ids/Diptera/Syrphidae.R')
make.labels(syrphid.list, spec=spec, table.name="syrphid.csv", genID="syrphid")

source('Final species ids/Lepidoptera/Moths.R')
make.labels(moths.morpho, spec=spec, table.name="moth.csv", genID="moth")

source('Species_ids/misc_flies.R')
make.labels(Misc_flies, spec=spec, table.name="missing2.csv", genID="bee")

