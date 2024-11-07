rm(list=ls())

## packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtreeExtra")
library(ggtreeExtra)
BiocManager::install("phyloseq")
library(phyloseq)

BiocManager::install("SparseArray")
library(SparseArray)



library(tidyr)
library(dplyr)
library(bipartite)
library(phyloseq)
library(TreeTools)
library(devtools)
library(ape)


BiocManager::install("TreeSummarizedExperiment")
library(TreeSummarizedExperiment)

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)


library(ggtree)
library(tidyverse)
library(ggtreeExtra)
library(ggstar)
library(ggplot2)
library(treeio)
library(ggnewscale)
library(tibble)
library(pals)
library(viridis)
library(phyloseq)
library(randomcoloR)
library(phytools)

## local directory based on person running
## local.path <- "C:/Users/rah10/"
local.path <- "~/"

wdpath <-
  file.path(local.path, 'Dropbox (University of Oregon)/skyIslands/data')

setwd(wdpath)

## Data imports
spec16s <- read.csv('spec_RBCL_16s.csv') %>%
  filter(Apidae==1) 

names(spec16s) <- gsub(x = names(spec16s), pattern = "X1",
                       replacement = "1")  
names(spec16s) <- gsub(x = names(spec16s), pattern = "RBCL.",
                       replacement = "RBCL:")  

meta_cols <- c('UniqueID', 'Family', 'Genus', 'Species',
               'GenusSpecies', 'Sex', 'Site', 'Meadow')

# meta <- spec16s %>%
#   select(all_of(meta_cols), Apidae, starts_with('16s')) %>%
#   na.omit() %>%
#   select(!starts_with('X16s')) %>%
#   filter(Apidae == 1) #%>%
# mutate(Site = factor(Site, levels=c("JC", ## ordered by latitude north to south
#                                     "SM",
#                                     "SC",
#                                     "MM",
#                                     "HM",
#                                     "PL",
#                                     "CH")))


## can probs figure out a way to do this automatically but was getting
## frustrated so hard coded this
apis_sites <- c('MM', 'HM', 'PL', 'CH')
bombus_sites <- c('JC', 'SM', 'SC', 'MM', 'HM', 'PL', 'CH')
anthophora_sites <- c('JC', 'SM', 'SC', 'MM', 'HM')
megachile_sites <- c('JC', 'SM', 'SC', 'MM', 'HM', 'PL')


## functions
##match unique IDs
match_shared_ID <- function(first_df, second_df){
  shared <- first_df$UniqueID[first_df$UniqueID %in% second_df$UniqueID]
  matched_df <- first_df %>%
    filter(UniqueID %in% shared)
  matched_df
}

#match tip labels
match_shared_tiplabels <- function(tree, pres_abs_table){ #input a phyloseq tree and a presence absense table
  tree_tips <- tree$tip.label #create object to store tree tip labels (strains)
  #browser()
  all_cols <- pres_abs_table %>% #filter the pres/abs table to remove unique ID
    select(-UniqueID) 
  #browser()
  match_cols <- all_cols[colnames(all_cols) %in% tree_tips] #match the features to the tree tips
  #browser()
  match_cols$UniqueID <- pres_abs_table$UniqueID # add back in the unique IDs
  match_cols #return the updated pres abs table with the tiplabels matched to the tree
}

#collapse tips that have the same label into one clade
collapse_identical_tips <- function(phy,tip_label){
  #matching_tips is initialized with the indices of tips in the phylogenetic tree (phy) whose labels match the provided tip_label. The function identifies all tips with the same label.
  matching_tips <- which(phy$tip.label==tip_label)
  nt <- length(phy$tip.label) # number of tips in tree
  nm <- length(matching_tips) # Number of tips matching the label
  keep <- numeric(nm) #keep is initialized as a numeric vector of length nm. It is used to keep track of which tips should be retained (1) and which tips should be dropped (0) in the new tree.
  
  #The while loop iterates through the indices of matching_tips to determine which tips to keep and which to drop.
  cur_tip <- 1
  #Inside the loop, the variable cur_tip is the current tip being considered, and next_tip is the tip immediately after cur_tip.
  while(cur_tip<=nm){
    if(cur_tip == nm){
      keep[cur_tip] <- 1
      break
    }
    next_tip <- cur_tip + 1
    #mrca_ (most recent common ancestor) is calculated using the getMRCA function for the tips identified by matching_tips[cur_tip] and matching_tips[next_tip]. This helps find the common ancestor of the current tip and the next tip.
    mrca_ <- getMRCA(phy,c(matching_tips[cur_tip],matching_tips[next_tip]))
    #descendants contains the indices of all descendants of the common ancestor, which includes both tips and internal nodes.
    descendants <- getDescendants(phy, mrca_)
    #descendant_tips is calculated to include only those indices from descendants that correspond to actual tips (i.e., indices less than or equal to nt).
    descendant_tips <- descendants[descendants<=nt]
    #The function checks if all descendant_tips are in the list of matching_tips. If they are, it means all these tips can be collapsed into a single branch, and they are marked to be kept.
    #The variable keep is updated accordingly, and cur_tip is advanced to skip the tips that have been collapsed.
    if(all(descendant_tips %in% matching_tips)){
      keep[cur_tip] <- 1
      cur_tip <- cur_tip + length(descendant_tips)
    }else{ #If not all descendant_tips are in the list of matching_tips, it means that not all tips can be collapsed, so the current tip is marked to be kept, and cur_tip is incremented by 1.
      keep[cur_tip] <- 1
      cur_tip <- cur_tip + 1
    }
  }
  #After the loop completes, the to_drop variable contains the indices of the tips that need to be dropped to collapse identical labels.
  to_drop <- matching_tips[!keep]
  #The function then creates a new phylogenetic tree (new_phy) by using the drop.tip function to remove the tips identified in to_drop.
  new_phy <- drop.tip(phy,to_drop)
  #browser()
  #Finally, the new phylogenetic tree is returned as the output of the function.
  return(new_phy)
}

##function that filters a tree object to a certain sample bee genus, generates a presence absence heatmap for whether that tree tip was found 
## at any given site, then plots the tree with the appended heatmap


##probs want to make sure each site in each individual graph is colored uniformly for all graphs 

phylotree_heatmap_byGenus <- function(tree.object, metadata, genus.or.spp, this.species, presAbsTable, site.order, all_levels=TRUE, levels_to_drop, clade_names=NULL, do_collapse=FALSE){
  #filter to include just the unique IDs in the specified genus
  if (genus.or.spp=='Species'){
    sp_ids <- metadata %>%
      filter(GenusSpecies==this.species) %>%
      select(UniqueID)
    
    #pull out all sites that include the specified genus
    my_sites <- unique(metadata$Site[metadata$GenusSpecies==this.species])
  }
  if (genus.or.spp=='Genus'){
    sp_ids <- metadata %>%
      filter(Genus==this.species) %>%
      select(UniqueID)
    
    #pull out all sites that include the specified genus
    my_sites <- unique(metadata$Site[metadata$Genus==this.species])
  }
  #browser()
  #remove tips from the tree that are not in the list of unique IDs in the specified genus
  trimmed_tree <- prune_samples(rownames(tree.object@sam_data) %in% sp_ids$UniqueID, tree.object)
  
  #remove taxa from the tree where there were less than zero observations
  pruned_tree <- prune_taxa(taxa_sums(trimmed_tree) > 0, trimmed_tree)
  
  ## for apis drop the incorrect bifidobacterium that is sorting with orbaceae
  #pruned_tree <- prune_taxa(!"16s:d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Bifidobacteriales; f__Bifidobacteriaceae; g__Bifidobacterium; s__Bifidobacterium_coryneforme", pruned_tree)
  
  
  #read in the taxonomic info for each feature
  feature.2.tax.16s <-
    read.table("SI_pipeline/merged/16s/taxonomy.tsv", sep="\t",
               header=TRUE)
  
  #make labels from the taxonomic info
  feature.2.tax.16s$Taxon <- paste("16s", feature.2.tax.16s$Taxon, sep=':')
  
  # ## convert to a phylo class which is more useful downstream
  gentree <- phy_tree(pruned_tree, errorIfNULL=TRUE)
  
  ## match the tip labs to the table with feature ID and Taxon
  gentree$tip.label  <-  feature.2.tax.16s$Taxon[match(gentree$tip.label,
                                                       feature.2.tax.16s$Feature.ID)]
  
  
  # Identify tips with labels exactly matching '16s:D_0__Bacteria'
  matching_tips <- grep('^16s:D_0__Bacteria$', gentree$tip.label)
  
  # Drop the matching tips
  gentree <- drop.tip(gentree, matching_tips)
  orig_len <- length(gentree$tip.label)
  print(paste('original # tips', orig_len))
  
  if (do_collapse == TRUE){
    #pull out unique tip labels
    groups <- unique(gentree$tip.label)
    
    #collapse branches that have the same label
    for (this.group in groups){
      gentree <- collapse_identical_tips(gentree, this.group)
    }
  }
  
  
  print(length(gentree$tip.label))
  if (do_collapse == TRUE){
    #pull out unique tip labels
    original_nodes <- gentree %>% as_tibble()
    groups <- unique(gentree$tip.label)
    
    
    # collapse branches that have the same label
    for (this.group in groups){
      gentree <- collapse_identical_tips(gentree, this.group)
    }
    new_nodes <- gentree %>% as_tibble()
  }
  print(length(gentree$tip.label))
  
  
  # if(all_levels==FALSE){
  #   if(final_level == ' s__'){
  #     rest_of_label <- 'to genus'
  #   } else if(final_level == ' g__'){
  #     rest_of_label <- 'to family'
  #   }
  #   this_level <- paste(": Collapsed", rest_of_label)
  # } else {this_level <- ': Full Tree'}
  
  
  matched_presabs <- match_shared_tiplabels(gentree, presAbsTable)
  ## dropping lots of tiplabels here...
  #browser()
  matched_pres_meta <- match_shared_ID(matched_presabs, metadata)
  
  
  
  matched_id <- matched_pres_meta$UniqueID
  row.names(matched_pres_meta) <- matched_id
  if (genus.or.spp=='Species'){
    meta_match_sites <- match_shared_ID(metadata, matched_pres_meta) %>%
      select(UniqueID, Site, GenusSpecies) %>%
      mutate(Site = factor(Site)) %>%
      filter(GenusSpecies==this.species) %>%
      select(!GenusSpecies) %>%
      group_by(UniqueID, Site) %>%
      mutate(n= n()) %>%
      pivot_wider(names_from=Site,
                  values_from = n,
                  names_expand = TRUE,
                  id_expand=TRUE) %>%
      pivot_longer(cols=2:length(colnames(.)),
                   names_to='Site',
                   values_to='Site_present') %>%
      filter(Site_present > 0) #%>%
    #mutate(Site = factor(Site, levels=site.order))
    #browser()
  }
  #browser()
  
  if (genus.or.spp=='Genus'){
    meta_match_sites <- match_shared_ID(metadata, matched_pres_meta) %>%
      select(UniqueID, Site, Genus) %>%
      mutate(Site = factor(Site)) %>%
      filter(Genus==this.species) %>%
      select(!Genus) %>%
      group_by(UniqueID, Site) %>%
      mutate(n= n()) %>%
      pivot_wider(names_from=Site,
                  values_from = n,
                  names_expand = TRUE,
                  id_expand=TRUE) %>%
      pivot_longer(cols=2:length(colnames(.)),
                   names_to='Site',
                   values_to='Site_present') %>%
      filter(Site_present > 0) #%>%
    #mutate(Site = factor(Site, levels=site.order))
    #browser()
  }
  
  
  features_site_metadata <- match_shared_ID(matched_pres_meta, meta_match_sites) %>%
    right_join(meta_match_sites, by='UniqueID') %>%
    #browser()
    pivot_longer(cols = starts_with('16s'), names_to = 'bacteria', values_to = 'bact_pres') %>%
    group_by(bacteria) %>%
    filter(bact_pres == 1) %>%
    select(!bact_pres) %>%
    relocate(bacteria) %>%
    group_by(bacteria) %>%
    add_count(Site, name="n_individuals") %>%
    mutate(SiteCount = as.numeric(n_distinct(Site))) %>%
    mutate(Obligate = as.numeric(str_detect(bacteria, "Lactobacillaceae|Bifidobacteriaceae|Neisseriaceae|Orbaceae|Bartonella|Acetobacteraceae")))
  #browser() 
  
  # Preprocess data to create a 'color' column
  features_site_metadata <- features_site_metadata %>%
    mutate(obligate_color = case_when(
      grepl("Lactobacillaceae", bacteria) ~ "#882D17",   # Lactobacillaceae
      grepl("Bifidobacteriaceae", bacteria) ~ "#B3446C", # Bifidobacteriaceae
      grepl("Neisseriaceae", bacteria) ~ "#DCD300",      # Neisseriaceae
      grepl("Orbaceae", bacteria) ~ "#8DB600",           # Orbaceae
      grepl("Bartonella", bacteria) ~ "#604E97",         # Bartonella
      grepl("Acetobacteraceae", bacteria) ~ "#F6A600",   # Acetobacteraceae
      TRUE ~ NA_character_  # Default NA for non-matching cases
    ))
  
  tree_tip_labs <- gentree$tip.label
  
  #dropping branches that weren't in the presence abs table
  final_drop <- gentree$tip.label[!(gentree$tip.label %in% features_site_metadata$bacteria)]
  #browser()
  
  if (length(final_drop) > 0){
    gentree <- drop.tip(gentree, final_drop)
  }
  
  ## save out order of tips
  is_tip <- gentree$edge[,2] <= length(gentree$tip.label)
  
  ordered_tips <- gentree$edge[is_tip, 2]
  
  tip.order <- gentree$tip.label[ordered_tips]
  
  
  p <- ggtree(gentree, layout='rectangular') 
  p
  
  p2 <- p +
    new_scale_fill() + 
    #geom_tiplab(align=TRUE, size=2) + 
    coord_cartesian(clip="off") +
    geom_fruit(
      data=features_site_metadata,
      geom=geom_tile,
      pwidth=0.05,
      offset=0.1,
      mapping=aes(y=bacteria,
                  #x=SiteCount,
                  fill=SiteCount, width=0.1),
      show.legend=TRUE) +
    labs(fill='Number of Sites')+
    scale_fill_gradient(high = "black", low ="lightgrey") +
    #ggtitle(paste(this.species, this_level)) +
    new_scale_fill() + 
    geom_fruit(
      data=features_site_metadata,
      geom=geom_tile,
      pwidth=0.05,
      offset=0.085,
      mapping=aes(y=bacteria,
                  #x=SiteCount,
                  fill=obligate_color, width=0.1),
      show.legend=FALSE) +
    scale_fill_identity() + # Use exact colors from the data
    new_scale_fill() + 
    #geom_tiplab(align=TRUE, linetype='dashed', aes(label = "")) +
    # geom_tippoint(aes(
    #  subset=(!grepl("Lactobacillus|Bifidobacterium|Snodgrassella|Gilliamella|Frischella|Bartonella|Commensalibacter",label,fixed=TRUE)==TRUE)), pch=15, color='black')+
    geom_tippoint(aes(
      subset=(grepl("Orbaceae",label,fixed=TRUE)==TRUE)), pch=21, fill="#8DB600", size=4)+
    geom_tippoint(aes(
      subset=(grepl("Lactobacillaceae",label,fixed=TRUE)==TRUE)), pch=21, fill="#882D17", size=4)+
    geom_tippoint(aes(
      subset=(grepl("Neisseriaceae",label,fixed=TRUE)==TRUE)), pch=21, fill="#DCD300", size=4)+
    geom_tippoint(aes(
      subset=(grepl( "Bifidobacteriaceae",label,fixed=TRUE)==TRUE)), pch=21, fill="#B3446C", size=4)+
    geom_tippoint(aes(
      subset=(grepl( "Acetobacteraceae",label,fixed=TRUE)==TRUE)), pch=21, fill="#F6A600", size=4)+
    geom_tippoint(aes(
      subset=(grepl("Bartonella",label,fixed=TRUE)==TRUE)), pch=21, fill="#604E97", size=4) 
  
  ## list [[1]] is tree, [[2]] is metadata, [[3]] is tip.order
  
  return(list(p2, features_site_metadata, tip.order))
  
  
}

