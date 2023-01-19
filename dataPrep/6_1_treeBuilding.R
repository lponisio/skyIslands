
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

spec16s <- read.csv('spec_RBCL_16s.csv')

meta_cols <- c('UniqueID', 'Family', 'Genus', 'Species', 'Sex', 'GeographyFK', 'Site', 'Meadow')

meta <- spec16s %>%
  select(all_of(meta_cols), Apidae) %>%
  filter(Apidae == 1)

### Need to filter phylogeny by uniqueID genus

# phylotree_by_genus <- function(tree.object, metadata, genus){
#   genus_ids <- metadata %>%
#     filter(Genus==genus) %>%
#     select(UniqueID)
#   
#   trimmed_tree <- prune_samples(sample_names(tree.object) %in% genus_ids$UniqueID, tree.object)
#   
#   ggtree(trimmed_tree, layout='rectangular')
# }
# 
# phylotree_by_genus(apis_phyloseq, meta, "Apis")
# phylotree_by_genus(bom)

apis_ids <- meta %>%
              filter(Genus=='Apis') %>%
              select(UniqueID)

bombus_ids <- meta %>%
  filter(Genus=='Bombus') %>%
  select(UniqueID)

anthophora_ids <- meta %>%
  filter(Genus=='Anthophora') %>%
  select(UniqueID)

megachile_ids <- meta %>%
  filter(Genus=='Megachile') %>%
  select(UniqueID)


####subsetting tree

apis_phyloseq <- prune_samples(sample_names(physeq16sR0) %in% apis_ids$UniqueID, physeq16sR0)
bombus_phyloseq <- prune_samples(sample_names(physeq16sR0) %in% bombus_ids$UniqueID, physeq16sR0)
anthophora_phyloseq <- prune_samples(sample_names(physeq16sR0) %in% anthophora_ids$UniqueID, physeq16sR0)
megachile_phyloseq <- prune_samples(sample_names(physeq16sR0) %in% megachile_ids$UniqueID, physeq16sR0)

###
#now the sample numbers are correct but it is still plotting the same tree for each genus
#probs want to use prune_taxa but need to first find out which taxa are still in subsetted genus uniqueID lists
#will need to access genus_ids and make a list of all remaining X16s files, then use similar code to how we selected 
#samples to prune filtered by genus
###

ggtree(anthophora_phyloseq, layout='rectangular')


# ##updating features to taxa
# setwd('../../skyIslands_saved')
# feature.2.tax.16s <-
#   read.table("SI_pipeline/merged/16s/taxonomy16s.txt", sep="\t",
#              header=TRUE)
# 
# feature.2.tax.16s$Taxon <- paste("16s", feature.2.tax.16s$Taxon, sep=':')
# 
# ## convert to a phylo class which is more useful downstream
# apis.tree.16sR0 <- phy_tree(apis_phyloseq, errorIfNULL=TRUE)
# 
# ## match the tip labs to the table with feature ID and Taxon
# apis.tree.16sR0$tip.label  <-  feature.2.tax.16s$Taxon[match(apis.tree.16sR0$tip.label,
#                                                         feature.2.tax.16s$Feature.ID)] 
# ##

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


comm_presabs <- as.data.frame(indiv.comm.16sR0) 
comm_presabs[comm_presabs > 0] <- 1 #converts from abundance to P/A -- check if needs to actualy do this


matched_presabs <- match_shared_tiplabels(apis.tree.16sR0, comm_presabs)
matched_pres_meta <- match_shared_ID(matched_presabs, meta)


# comm_presabs <- as.data.frame(indiv.comm.16sR0) 
# comm_presabs[comm_presabs > 0] <- 1 #converts from abundance to P/A -- check if needs to actualy do this

comm_presabs <- tibble::rownames_to_column(comm_presabs, "UniqueID")
matched_pres_meta <- match_shared_ID(comm_presabs, meta)
matched_id <- matched_pres_meta$UniqueID
row.names(matched_pres_meta) <- matched_id

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
  mutate(Site = factor(Site, levels=c("JC", ## ordered by latitude north to south
                                      "SM",
                                      "SC",
                                      "MM",
                                      "HM",
                                      "PL",
                                      "CH")))



features_genus_metadata <- match_shared_ID(matched_pres_meta, meta_match_genus) %>%
  right_join(meta_match_genus, by='UniqueID') %>%
  pivot_longer(cols = starts_with('16s'), names_to = 'bacteria', values_to = 'bact_pres') %>%
  group_by(bacteria) %>%
  filter(bact_pres == 1) %>%
  select(!bact_pres) %>%
  relocate(bacteria) %>% 
  mutate(Genus = factor(Genus, levels=c("Apis", ## ordered by relative (eyeballed) PCA microbiome similarity
                                      "Bombus",
                                      "Anthophora",
                                      "Megachile",
                                      "Agapostemon")))

bar_metadata <- features_genus_metadata %>%
  group_by(bacteria, Genus) %>%
  summarize(Count = sum(Genus_present)) %>%
  filter(Count>0)

## probs want to use trimmed tree for final
## visualization but right now will retain


p <- ggtree(tree.16s, layout='rectangular') 
p

# p <- ggtree(tree.16s,
#             layout="fan",
#             open.angle=10,
#             size=0.5,
#             aes(color=features_genus_metadata$Genus[]))
# 

p2 <- p + 
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


