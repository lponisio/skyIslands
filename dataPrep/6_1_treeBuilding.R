

##making phytools phylogeny tutorial http://www.phytools.org/Cordoba2017/ex/2/Intro-to-phylogenies.html

# library(ape)
# library(phangorn)
# library(phytools)
# library(geiger)
# library(devtools)
# install_github("liamrevell/phytools")
# 
# 
# taxonomy.table <- physeq16sR0@tax_table
# 
# plotTree(tree.16sR0,type="fan",fsize=0.7,lwd=1,
#          ftype="off") ##not great

###tutorial https://yulab-smu.top/treedata-book/chapter4.html



library(treeio)
library(ggtree)

# tree.16s$tip.label <- gsub("16s:", "", tree.16s$tip.label)
# bees.16s
# 
# ggtree(tree.16sR0, branch.length='none', layout='circular')
# 
# groupInfo <- data.frame(physeq16sR0@tax_table)

ggtree(tree.16s, layout='circular')


## import metadata

wdpath <- 'C:/Users/rah10/Dropbox (University of Oregon)/PonisioLocal/skyIslands'
setwd(wdpath)

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

# install.packages("remotes")
# remotes::install_github("xiangpin/MicrobiotaProcess")
# 
# #apis
# apis_subset <- meta %>%
#   filter(Genus == 'Apis') 
# apis_ids <- apis_subset %>%
#   select(UniqueID) 
# apis_taxonomy <- indiv.comm.16sR0[row.names(indiv.comm.16sR0) %in% apis_ids$UniqueID == TRUE,] 
# 
# apis_phylo <- as.phyloseq(apis_taxonomy)
# 
# ggtree(as.data.frame(apis_taxonomy), branch.length='none', layout='circular')
# 
# #bombus
# bombus_subset <- meta %>%
#   filter(Genus == 'Bombus') 
# bombus_ids <- bombus_subset %>%
#   select(UniqueID) 
# bombus_taxonomy <- indiv.comm.16sR0[row.names(indiv.comm.16sR0) %in% bombus_ids$UniqueID == TRUE,]
# 
# #megachile
# megachile_subset <- meta %>%
#   filter(Genus == 'Megachile') 
# megachile_ids <- megachile_subset %>%
#   select(UniqueID) 
# megachile_taxonomy <- indiv.comm.16sR0[row.names(indiv.comm.16sR0) %in% megachile_ids$UniqueID == TRUE,] 
# 
# #anthophora
# anthophora_subset <- meta %>%
#   filter(Genus == 'Anthophora') 
# anthophora_ids <- anthophora_subset %>%
#   select(UniqueID) 
# anthophora_taxonomy <- indiv.comm.16sR0[row.names(indiv.comm.16sR0) %in% anthophora_ids$UniqueID == TRUE,] 

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


spec16s <- read.csv('spec_RBCL_16s.csv')

meta_cols <- c('UniqueID', 'Family', 'Genus', 'Species', 'Sex', 'GeographyFK', 'Site', 'Meadow')

microbes <- spec16s %>%
  select(UniqueID, Site, Genus, starts_with('X16s')) %>%
  filter(!Genus == 'Agapostemon') %>%
  na.omit()

meta <- spec16s %>%
  select(all_of(meta_cols), Apidae) %>%
  filter(Apidae == 1)


##########################

##match unique IDs
match_shared_ID <- function(first_df, second_df){
  shared <- first_df$UniqueID[first_df$UniqueID %in% second_df$UniqueID]
  matched_df <- first_df %>%
  filter(UniqueID %in% shared)
  matched_df
}

match_shared_tiplabels <- function(tree, pres_abs_table){
  tree_tips <- tree$tip.label
  #browser()
  all_cols <- pres_abs_table %>%
    select(-UniqueID) 
  #browser()
  match_cols <- all_cols[colnames(all_cols) %in% tree_tips]
  #browser()
  match_cols$UniqueID <- pres_abs_table$UniqueID
  match_cols
}

matched_presabs <- match_shared_tiplabels(tree.16s, comm_presabs)
matched_pres_meta <- match_shared_ID(matched_presabs, meta)

#hmm rn there are more tip.labels in the tree than feature columns
# in the metadata... why???
#want to compare tree.16s$tip.label(265) to colnames in comm_presabs(254)

#trying to find the step that makes the numbers not match
#tree.16sR0
p <- ggtree(tree.16sR0, layout='circular')
n_occur <- data.frame(table(tree.16sR0$tip.label))

#253 features but some features have multiple occurences
# need to figure out a way to pair up the unique features
# with metadata -- ideally i want rows as tip.names
# with each column as metadata
# will need a column of presence and absence for
# each site
# each genus
comm_presabs <- as.data.frame(indiv.comm.16sR0) 
comm_presabs[comm_presabs > 0] <- 1 #converts from abundance to P/A
library(tibble)
comm_presabs <- tibble::rownames_to_column(comm_presabs, "UniqueID")
matched_pres_meta <- match_shared_ID(comm_presabs, meta)
matched_id <- matched_pres_meta$UniqueID
row.names(matched_pres_meta) <- matched_id
#matched_matrix <- matched_pres_meta %>% 
#  select(-UniqueID) %>%
#  t()

#making presence/abs columns in metadata
meta_match_sites <- match_shared_ID(meta, matched_pres_meta) %>%
  select(UniqueID, Site) %>%
  mutate(Site = as.factor(Site)) %>%
  group_by(UniqueID, Site) %>%
  count() %>%
  pivot_wider(UniqueID, 
              names_from=Site, 
              values_from = n, 
              values_fill=0,
              names_expand = TRUE,
              id_expand=TRUE) %>%
  pivot_longer(cols=CH:SM,
               names_to='Site',
               values_to='Site_present')


meta_match_genus <- match_shared_ID(meta, matched_pres_meta) %>%
  select(UniqueID, Genus) %>%
  mutate(Site = as.factor(Genus)) %>%
  group_by(UniqueID, Genus) %>%
  count() %>%
  pivot_wider(UniqueID, 
              names_from=Genus, 
              values_from = n, 
              values_fill=0,
              names_expand = TRUE,
              id_expand=TRUE) %>%
  pivot_longer(cols=Agapostemon:Megachile,
               names_to='Genus',
               values_to='Genus_present')

## now need to figure out how to incorporate
## the tip labels so that the features are 
## tied to the metadata -- features as rows
## columns for metadata

##matched_pres_meta is a table with unique ID 
## as rows and features as columns while both meta_match
## have unique ID as rows and metadata as columns
## need to create an item that has features as rows
## and metadata as columns



features_site_metadata <- match_shared_ID(matched_pres_meta, meta_match_sites) %>%
  right_join(meta_match_sites, by='UniqueID') %>%
  pivot_longer(cols = starts_with('16s'), names_to = 'bacteria', values_to = 'bact_pres') %>%
  group_by(bacteria) %>%
  filter(bact_pres == 1) %>%
  select(!bact_pres) %>%
  relocate(bacteria) %>% 
  mutate(strain_bee = paste(bacteria,UniqueID,sep="_"))

features_genus_metadata <- match_shared_ID(matched_pres_meta, meta_match_genus) %>%
  right_join(meta_match_genus, by='UniqueID') %>%
  pivot_longer(cols = starts_with('16s'), names_to = 'bacteria', values_to = 'bact_pres') %>%
  group_by(bacteria) %>%
  filter(bact_pres == 1) %>%
  select(!bact_pres) %>%
  relocate(bacteria) %>% 
  mutate(strain_bee = paste(bacteria,UniqueID,sep="_"))

## probs want to use trimmed tree for final
## visualization but right now will retain
## duplicates
#library(wesanderson)


p_rectangular <- ggtree(tree.16s, #layout="fan", open.angle=10, size=0.5)
            layout='rectangular')

p_fan <- ggtree(tree.16s, layout="fan", open.angle=10, size=0.5)

                p1

p2 <- p + 
  geom_fruit(
    data=features_site_metadata,
    geom=geom_tile,
    mapping=aes(y=bacteria, 
                x=Site, 
                alpha=as.factor(Site_present),
                fill=Site)) +
  scale_fill_manual(
    values=c("#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7"))
p2

p3 <- p + 
  geom_fruit(
    data=features_genus_metadata,
    geom=geom_tile,
    mapping=aes(y=bacteria, 
                x=Genus, 
                alpha=as.factor(Genus_present),
                fill=Genus)) 
p3