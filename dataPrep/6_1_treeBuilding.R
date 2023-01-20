
## tutorial: https://bioconductor.org/packages/devel/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html

#run up to line 61 in 6_rbcl_16sPrep.R

## packages
library(ggtree)
library(tidyverse)
library(ggtreeExtra)
library(ggstar)
library(ggplot2)
library(treeio)
library(ggnewscale)
library(tibble)
library(pals)

## working dir

wdpath <- 'C:/Users/rah10/Dropbox (University of Oregon)/PonisioLocal/skyIslands'
setwd(wdpath)

## Data imports

spec16s <- read.csv('spec_RBCL_16s.csv') %>%
  filter(Apidae==1) 

names(spec16s) <- gsub(x = names(spec16s), pattern = "X1", replacement = "1")  

meta_cols <- c('UniqueID', 'Family', 'Genus', 'Species', 'Sex', 'GeographyFK', 'Site', 'Meadow')

meta <- spec16s %>%
  select(all_of(meta_cols), Apidae, starts_with('16s')) %>%
  na.omit() %>%
  select(!starts_with('X16s')) %>%
  filter(Apidae == 1)

### Need to filter phylogeny by uniqueID genus
setwd("../../skyIslands_saved")

phylotree_by_genus <- function(tree.object, metadata, genus){
  genus_ids <- metadata %>%
    filter(Genus==genus) %>%
    select(UniqueID)

  trimmed_tree <- prune_samples(sample_names(tree.object) %in% genus_ids$UniqueID, tree.object)
  
  pruned_tree <- prune_taxa(taxa_sums(trimmed_tree) > 0, trimmed_tree)
  
  
  feature.2.tax.16s <-
    read.table("SI_pipeline/merged/16s/taxonomy16s.txt", sep="\t",
               header=TRUE)
  
  feature.2.tax.16s$Taxon <- paste("16s", feature.2.tax.16s$Taxon, sep=':')
  
  # ## convert to a phylo class which is more useful downstream
  gentree <- phy_tree(pruned_tree, errorIfNULL=TRUE)
  
  ## match the tip labs to the table with feature ID and Taxon
  gentree$tip.label  <-  feature.2.tax.16s$Taxon[match(gentree$tip.label,
                                                        feature.2.tax.16s$Feature.ID)]

  p <- ggtree(gentree, layout='rectangular')
  p
}

apis_tree <- phylotree_by_genus(physeq16sR0, meta, "Apis")
bombus_tree <- phylotree_by_genus(physeq16sR0, meta, "Bombus")
anthophora_tree <- phylotree_by_genus(physeq16sR0, meta, "Anthophora")
megachile_tree <- phylotree_by_genus(physeq16sR0, meta, "Megachile")


##########################

##match unique IDs
match_shared_ID <- function(first_df, second_df){
  shared <- first_df$UniqueID[first_df$UniqueID %in% second_df$UniqueID]
  matched_df <- first_df %>%
  filter(UniqueID %in% shared)
  matched_df
}
#######################################
#match tip labels
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

####
#community presence absence df

comm_presabs <- as.data.frame(indiv.comm.16sR0) 
comm_presabs[comm_presabs > 0] <- 1 
comm_presabs <- tibble::rownames_to_column(comm_presabs, "UniqueID")


#
make_heatmap_tree <- function(tree, presabs_table, metadata){
  
  matched_presabs <- match_shared_tiplabels(tree, presabs_table)
  
  matched_pres_meta <- match_shared_ID(matched_presabs, metadata)
  
  matched_id <- matched_pres_meta$UniqueID
  
  row.names(matched_pres_meta) <- matched_id
  
  #making presence/abs columns in metadata
  meta_match_sites <- match_shared_ID(metadata, matched_pres_meta) %>%
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
  
  features_site_metadata <- match_shared_ID(matched_pres_meta, meta_match_sites) %>%
    right_join(meta_match_sites, by='UniqueID') %>%
    pivot_longer(cols = starts_with('16s'), names_to = 'bacteria', values_to = 'bact_pres') %>%
    group_by(bacteria) %>%
    filter(bact_pres == 1) %>%
    select(!bact_pres) %>%
    relocate(bacteria) %>% 
    mutate(Site = factor(Site, levels=c("JC", ## ordered by latitude north to south
                                        "SM",
                                        "SC",
                                        "MM",
                                        "HM",
                                        "PL",
                                        "CH")))
  browser()
  
  p2 <- tree + 
    geom_fruit(
      data=features_site_metadata,
      geom=geom_tile,
      mapping=aes(y=bacteria, 
                  x=Site, 
                  alpha=Site_present,
                  fill=Site), 
      axis.params=list(
        axis="x",
        title = "Site",
        text.size=1.5,
        vjust=-99,
        #text.angle=-45
      ),
      show.legend=FALSE) 
  p2
  
  
  
}


apis_heat_tree <- make_heatmap_tree(apis_tree, comm_presabs, meta)
apis_heat_tree


#
comm_presabs <- tibble::rownames_to_column(comm_presabs, "UniqueID")
matched_pres_meta <- match_shared_ID(comm_presabs, meta)
matched_id <- matched_pres_meta$UniqueID
row.names(matched_pres_meta) <- matched_id
#
# #making presence/abs columns in metadata
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
#
#
# meta_match_genus <- match_shared_ID(meta, matched_pres_meta) %>%
#   select(UniqueID, Genus) %>%
#   mutate(Site = as.factor(Genus)) %>%
#   group_by(UniqueID, Genus) %>%
#   count() %>%
#   pivot_wider(UniqueID,
#               names_from=Genus,
#               values_from = n,
#               values_fill=0,
#               names_expand = TRUE,
#               id_expand=TRUE) %>%
#   pivot_longer(cols=Agapostemon:Megachile,
#                names_to='Genus',
#                values_to='Genus_present')

## now need to figure out how to incorporate
## the tip labels so that the features are
## tied to the metadata -- features as rows
## columns for metadata

##matched_pres_meta is a table with unique ID
## as rows and features as columns while both meta_match
## have unique ID as rows and metadata as columns
## need to create an item that has features as rows
## and metadata as columns


#
features_site_metadata <- match_shared_ID(matched_pres_meta, meta_match_sites) %>%
  right_join(meta_match_sites, by='UniqueID') %>%
  pivot_longer(cols = starts_with('16s'), names_to = 'bacteria', values_to = 'bact_pres') %>%
  group_by(bacteria) %>%
  filter(bact_pres == 1) %>%
  select(!bact_pres) %>%
  relocate(bacteria) %>%
  mutate(Site = factor(Site, levels=c("JC", ## ordered by latitude north to south
                                      "SM",
                                      "SC",
                                      "MM",
                                      "HM",
                                      "PL",
                                      "CH")))
#
#
#
# features_genus_metadata <- match_shared_ID(matched_pres_meta, meta_match_genus) %>%
#   right_join(meta_match_genus, by='UniqueID') %>%
#   pivot_longer(cols = starts_with('16s'), names_to = 'bacteria', values_to = 'bact_pres') %>%
#   group_by(bacteria) %>%
#   filter(bact_pres == 1) %>%
#   select(!bact_pres) %>%
#   relocate(bacteria) %>%
#   mutate(Genus = factor(Genus, levels=c("Apis", ## ordered by relative (eyeballed) PCA microbiome similarity
#                                       "Bombus",
#                                       "Anthophora",
#                                       "Megachile",
#                                       "Agapostemon")))
#
# bar_metadata <- features_genus_metadata %>%
#   group_by(bacteria, Genus) %>%
#   summarize(Count = sum(Genus_present)) %>%
#   filter(Count>0)

## probs want to use trimmed tree for final
## visualization but right now will retain


p <- ggtree(apistree, layout='rectangular', tip.label=FALSE) 
p


# 

p2 <- anthophora_tree + 
  geom_fruit(
    data=features_site_metadata,
    geom=geom_tile,
    mapping=aes(y=bacteria, 
                x=Site, 
                alpha=Site_present,
                fill=Site), 
    axis.params=list(
      axis="x",
      title = "Site",
      text.size=1.5,
      vjust=-99,
      #text.angle=-45
      )) 
p2


p3 <- p2 +
  #new_scale_fill() +
  geom_fruit(data = features_genus_metadata, 
             geom = geom_tile, 
             mapping=aes(y=bacteria,
                         x= Genus, 
                         alpha = Genus_present, 
                         #fill = Genus
                         ), 
             axis.params=list(
               axis="x",
               title = "Bee Genus",
               text.size=1,
               #text.angle=45,
               vjust=-125,
               #hjust=-10
           ), show.legend=FALSE) 

p3  


##heatmap version
p4 <- p2 + new_scale_fill()
gheatmap(p4, df2, offset=15, width=.3,
         colnames_angle=90, colnames_offset_y = .25) +
  scale_fill_viridis_c(option="A", name="continuous\nvalue")


