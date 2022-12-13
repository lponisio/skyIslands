

##making phytools phylogeny tutorial http://www.phytools.org/Cordoba2017/ex/2/Intro-to-phylogenies.html

library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(devtools)
install_github("liamrevell/phytools")

plotTree(tree.16sR0,type="fan",fsize=0.7,lwd=1,
         ftype="off") ##not great

###tutorial https://yulab-smu.top/treedata-book/chapter4.html



library(treeio)
library(ggtree)

ggtree(tree.16sR0, aes(color=meta$Site), branch.length='none', layout='circular')

## import metadata

spec16s <- read.csv('spec_RBCL_16s.csv')

meta_cols <- c('UniqueID', 'Family', 'Genus', 'Species', 'Sex', 'GeographyFK', 'Site', 'Meadow')

microbes <- spec16s %>%
  select(UniqueID, Site, Genus, starts_with('X16s')) %>%
  filter(!Genus == 'Agapostemon') %>%
  na.omit()

meta <- spec16s %>%
  select(all_of(meta_cols), Apidae) %>%
  filter(Apidae == 1)

#######################################################
##create lists of unique IDs for each species to filter
######################################################

#apis
apis_subset <- meta %>%
  filter(Genus == 'Apis') 
apis_ids <- apis_subset %>%
  select(UniqueID) 
apis_taxonomy <- taxonomy16sR0 %>%
  t() 
apis_taxonomy <- apis_taxonomy[row.names(apis_taxonomy) %in% apis_ids$UniqueID == TRUE,] %>%
  t()
  

#bombus
bombus_subset <- meta %>%
  filter(Genus == 'Bombus') 
bombus_ids <- bombus_subset %>%
  select(UniqueID) 
bombus_taxonomy <- taxonomy16sR0 %>%
  t() 
bombus_taxonomy <- bombus_taxonomy[row.names(bombus_taxonomy) %in% bombus_ids$UniqueID == TRUE,] %>%
  t()

#megachile
megachile_subset <- meta %>%
  filter(Genus == 'Megachile') 
megachile_ids <- megachile_subset %>%
  select(UniqueID) 
megachile_taxonomy <- taxonomy16sR0 %>%
  t() 
megachile_taxonomy <- megachile_taxonomy[row.names(megachile_taxonomy) %in% megachile_ids$UniqueID == TRUE,] %>%
  t()

#anthophora
anthophora_subset <- meta %>%
  filter(Genus == 'Anthophora') 
anthophora_ids <- anthophora_subset %>%
  select(UniqueID) 
anthophora_taxonomy <- taxonomy16sR0 %>%
  t() 
anthophora_taxonomy <- anthophora_taxonomy[row.names(anthophora_taxonomy) %in% anthophora_ids$UniqueID == TRUE,] %>%
  t()
