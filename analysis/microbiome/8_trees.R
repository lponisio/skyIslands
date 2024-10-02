
## tutorial: https://bioconductor.org/packages/devel/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html


rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path
setwd(local.path)
#####
setwd('skyIslands/')


source('dataPrep/src/trees_init.R')
load("indiv.comm16sR0.Rdata")
load('presAbsTable.Rdata')
load('spec_RBCL_16s.Rdata')
load('physeq16s.Rdata')
load('networks/microNets.Rdata')
load('spec_RBCL_16s.Rdata')
##import community presence/absence file
setwd("../../skyIslands_saved")

library(pals)
library(cowplot)
library(ggpubr)
library(gridExtra)


#######################################


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



meta_cols <- c('UniqueID', 'Family', 'Genus', 'Species', 'GenusSpecies', 'Sex', 'Site', 'Meadow')

meta <- spec.net %>%
  filter(Apidae == 1) %>%
    select(all_of(meta_cols), Apidae, starts_with('16s')) %>%
    na.omit() %>%
    select(!starts_with('16s')) 

comm_presabs <- as.data.frame(indiv.comm.16sR0) #load in the pres/abs table
comm_presabs[comm_presabs > 0] <- 1 #change all rel abund to 1
comm_presabs <- tibble::rownames_to_column(comm_presabs, "UniqueID") #make rownames (UniqueID) into column

finalASV <- as.data.frame(finalASVtable)
finalASV[finalASV > 0] <- 1 #change all rel abund to 1
finalASV <- tibble::rownames_to_column(finalASV, "UniqueID") #make rownames (UniqueID) into column


## Genus trees
# apis_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Apis", genus.or.spp='Genus', finalASV, apis_sites, all_levels=TRUE, do_collapse = TRUE)
# panelB <- apis_tree[[1]] + labs(tag="B. Apis")
# apis_meta <- apis_tree[[2]]
# panelB



# melissodes 
melissodes_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Melissodes", genus.or.spp='Genus', finalASV, melissodes_sites, do_collapse = TRUE)
panelB <- melissodes_tree[[1]] + labs(tag="B. Melissodes")
melissodes_meta <- melissodes_tree[[2]]
panelB



#bombus tree
bombus_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus", genus.or.spp='Genus', finalASV, bombus_sites, do_collapse = TRUE)
panelA <- bombus_tree[[1]] + labs(tag="A. Bombus")
bombus_meta <- bombus_tree[[2]]
panelA




## make legend


# Load Packages 
library("ggplot2") 
library("grid") 
library("gridExtra") 
library("cowplot") 

## obligate symbionts
these_obligates <- c("Acetobacteraceae",
                     "Bartonellaceae",
                     "Bifidobacteriaceae",
                     "Lactobacillaceae",
                     "Neisseriaceae",
                     "Orbaceae")


leg_col_order <- c("#F6A600", #Acetobacteriaceae
                   "#604E97", #Bartonellaceae
                   "#B3446C", #Bifidobacteriaceae
                   "#882D17", #Lactobacillaceae
                   "#DCD300", #Neisseriaceae
                   "#8DB600") #Orbaceae



# Create a DataFrame 
data.leg <- data.frame( 
  Xdata = rnorm(6),
  Ydata = rnorm(6), 
  Family = these_obligates,
  leg_color = leg_col_order)

# Create a Scatter Plot 
gplot <- ggplot(data.leg, aes(Xdata, Ydata, color = Family)) +    
  geom_point(size = 7) +
  scale_color_manual(values=data.leg$leg_color) +
  theme(legend.position='bottom') +
  labs(color='Bacteria Family') +
  guides(colour = guide_legend(nrow = 1)) + theme(legend.key=element_blank(),
                                                   legend.text=element_text(size=12),
                                                  legend.title=element_text(size=12))

# Draw Only Legend without plot 
# Grab legend from gplot 
panelD  <- get_legend(gplot)                     
plot(get_legend(gplot) )


#Set up the layout matrix so that the bottom legend spans both columns
layout <- rbind(c(1, 2),
                c(3, 3)) # The legend will span both columns

# Create the final layout
final_plot <- arrangeGrob(panelA, panelB, panelD,
                          layout_matrix = layout,
                          heights = c(9, 1)) # Adjust the heights as needed

# Center the legend panelD properly within its grid
panelD_centered <- arrangeGrob(panelD, 
                               ncol = 1, 
                               padding = unit(1, "lines"))  # Add padding if necessary

# Open a PDF device to save the plot
pdf("../skyIslands/analysis/microbiome/figures/final/grid_trees.pdf",
    height=8, width=11)  # Adjust width and height as needed

# Plot the final combined figure
grid.arrange(final_plot)

# Close the PDF device to complete saving
dev.off()



