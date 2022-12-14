

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

install.packages("remotes")
remotes::install_github("xiangpin/MicrobiotaProcess")

#apis
apis_subset <- meta %>%
  filter(Genus == 'Apis') 
apis_ids <- apis_subset %>%
  select(UniqueID) 
apis_taxonomy <- indiv.comm.16sR0[row.names(indiv.comm.16sR0) %in% apis_ids$UniqueID == TRUE,] 

apis_phylo <- as.phyloseq(apis_taxonomy)

ggtree(as.data.frame(apis_taxonomy), branch.length='none', layout='circular')

#bombus
bombus_subset <- meta %>%
  filter(Genus == 'Bombus') 
bombus_ids <- bombus_subset %>%
  select(UniqueID) 
bombus_taxonomy <- indiv.comm.16sR0[row.names(indiv.comm.16sR0) %in% bombus_ids$UniqueID == TRUE,]

#megachile
megachile_subset <- meta %>%
  filter(Genus == 'Megachile') 
megachile_ids <- megachile_subset %>%
  select(UniqueID) 
megachile_taxonomy <- indiv.comm.16sR0[row.names(indiv.comm.16sR0) %in% megachile_ids$UniqueID == TRUE,] 

#anthophora
anthophora_subset <- meta %>%
  filter(Genus == 'Anthophora') 
anthophora_ids <- anthophora_subset %>%
  select(UniqueID) 
anthophora_taxonomy <- indiv.comm.16sR0[row.names(indiv.comm.16sR0) %in% anthophora_ids$UniqueID == TRUE,] 
