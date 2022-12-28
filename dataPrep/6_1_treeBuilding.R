

##making phytools phylogeny tutorial http://www.phytools.org/Cordoba2017/ex/2/Intro-to-phylogenies.html

library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(devtools)
install_github("liamrevell/phytools")


taxonomy.table <- physeq16sR0@tax_table

plotTree(tree.16sR0,type="fan",fsize=0.7,lwd=1,
         ftype="off") ##not great

###tutorial https://yulab-smu.top/treedata-book/chapter4.html



library(treeio)
library(ggtree)

tree.16s$tip.label <- gsub("16s:", "", tree.16s$tip.label)
bees.16s

ggtree(tree.16sR0, branch.length='none', layout='circular')

groupInfo <- data.frame(physeq16sR0@tax_table)

ggtree(tree.16s, layout='circular', branch.length = 'none')

data(chiroptera, package="ape")
groupInfo <- split(chiroptera$tip.label, gsub("_///w+", "", chiroptera$tip.label))
chiroptera <- groupOTU(chiroptera, groupInfo)
ggtree(chiroptera, aes(color=group), layout='circular') + geom_tiplab(size=1, aes(angle=angle))


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

################## 12/27/22
## tutorial: https://bioconductor.org/packages/devel/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html

library(ggtreeExtra)
library(ggstar)
library(ggplot2)
library(ggtree)
library(treeio)
library(ggnewscale)

## what we need to do:
# 1. make tree object from 6_rbcl_16sPrep.R
# 2. get metadata
## use bees.16s and spec_16s_rbcl or whatever
####a. we want site, genus, species(?)
## tips in indiv.comm.16sR0
## make a table of presence and absence in indiv.comm.16sR0 for each unique ID
## then for each metadata (site and bee genus) attach those to the unique ID
## the dimensions are not the same, so will have to account somehow for 
## duplicate features (strains) so that the metadata dimensions are the 
## same and # of tree.16s tips

