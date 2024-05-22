
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
  #remove tips from the tree that are not in the list of unique IDs in the specified genus
  trimmed_tree <- prune_samples(sample_names(tree.object) %in% sp_ids$UniqueID, tree.object)
  
  #remove taxa from the tree where there were less than zero observations
  pruned_tree <- prune_taxa(taxa_sums(trimmed_tree) > 0, trimmed_tree)
  
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
  
  
  
# if(all_levels==FALSE){
# for (level in levels_to_drop) {
#       tips_to_drop <- grep(level, gentree$tip.label)
#       print(paste('tip drop len', length(tips_to_drop)))
#       if (length(tips_to_drop) > 0) {
#         gentree <- drop.tip(gentree, tips_to_drop)
#         new_len <- length(gentree$tip.label)
#         percent_drop <- (length(tips_to_drop)/orig_len)*100
#         print(paste(level, "percent tips dropped", percent_drop))
#         final_level <- level
#       }
#     }
# }
  if (all_levels == FALSE) {
    for (level in levels_to_drop) {
      #make a list of tip labels that include the specified level
      tips_to_modify <- grep(level, gentree$tip.label)
      orig.tip.location <- as_tibble(gentree) %>%
        select(node, label) %>%
        mutate(do_modify = {ifelse(.$node %in% tips_to_modify, 1, 0)})
      #browser()
      if (length(tips_to_modify) > 0) {
        # Replace the pattern and everything following it with an empty string
        gentree$tip.label[tips_to_modify] <- gsub(paste0(";", level, ".*"), "", gentree$tip.label[tips_to_modify])
        this.tip.level.label <- paste(level)
        new.tip.location <- orig.tip.location %>%
          mutate(!!this.tip.level.label := {ifelse(do_modify == 1, gsub(paste0(";", level, ".*"), "", label), label)})
        new_len <- length(gentree$tip.label)
        percent_drop <- (length(tips_to_modify) / orig_len) * 100
        print(paste(level, "percent tips dropped", percent_drop))
        final_level <- level
      }
      }
  }
  
  #browser()
    #print(length(gentree$tip.label))
  if (do_collapse == TRUE){
    #pull out unique tip labels
    original_nodes <- gentree %>% as_tibble() 
    groups <- unique(gentree$tip.label)
    #browser()
    
    #collapse branches that have the same label
    for (this.group in groups){
      gentree <- collapse_identical_tips(gentree, this.group)
    }
    new_nodes <- gentree %>% as_tibble()
  }
  #print(length(gentree$tip.label))
  #browser()

  if(all_levels==FALSE){
    if(final_level == ' s__'){
      rest_of_label <- 'to genus'
    } else if(final_level == ' g__'){
      rest_of_label <- 'to family'
    }
    this_level <- paste(": Collapsed", rest_of_label)
  } else {this_level <- ': Full Tree'}

  
  matched_presabs <- match_shared_tiplabels(gentree, presAbsTable)
  
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
    count() %>%
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
  if (genus.or.spp=='Genus'){
    meta_match_sites <- match_shared_ID(metadata, matched_pres_meta) %>%
      select(UniqueID, Site, Genus) %>%
      mutate(Site = factor(Site)) %>%
      filter(Genus==this.species) %>%
      select(!Genus) %>%
      group_by(UniqueID, Site) %>%
      count() %>%
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
    mutate(Obligate = as.numeric(str_detect(bacteria, "Lactobacillus|Bifidobacterium|Snodgrassella|Gilliamella|Frischella|Bartonella|Commensalibacter")))
   #browser() 
  tree_tip_labs <- gentree$tip.label
  
  #dropping branches that weren't in the presence abs table
  final_drop <- gentree$tip.label[!(gentree$tip.label %in% features_site_metadata$bacteria)]
#browser()
  
  if (length(final_drop) > 0){
    gentree <- drop.tip(gentree, final_drop)
  }
  
  p <- ggtree(gentree, layout='rectangular') 
  p
  
  p2 <- p +
    new_scale_fill() + 
    #geom_tiplab(align=TRUE, size=2) + 
    coord_cartesian(clip="off") +
    geom_fruit(
      data=features_site_metadata,
      geom=geom_tile,
      pwidth=0.1,
      offset=0.2,
      mapping=aes(y=bacteria,
                  #x=SiteCount,
                  fill=SiteCount, width=0.1),
      show.legend=TRUE) +
    labs(fill='Number of Sites')+
    scale_fill_gradient(high = "black", low ="lightgrey") +
    #ggtitle(paste(this.species, this_level)) +
    new_scale_fill() + 
    #geom_tippoint(aes(
    #  subset=(!grepl("Lactobacillus|Bifidobacterium|Snodgrassella|Gilliamella|Frischella|Bartonella|Commensalibacter",label,fixed=TRUE)==TRUE)), pch=15, color='black')+
    geom_tippoint(aes(
      subset=(grepl("Gilliamella",label,fixed=TRUE)==TRUE)), pch=21, fill="#E69F00", size=3)+
    geom_tippoint(aes(
      subset=(grepl("Lactobacillus",label,fixed=TRUE)==TRUE)), pch=21, fill="#56B4E9", size=3)+
    geom_tippoint(aes(
      subset=(grepl("Snodgrassella",label,fixed=TRUE)==TRUE)), pch=21, fill="#009E73", size=3)+
    geom_tippoint(aes(
      subset=(grepl("Commensalibacter",label,fixed=TRUE)==TRUE)), pch=21, fill="#F0E442", size=3)+
    geom_tippoint(aes(
      subset=(grepl( "Bifidobacterium",label,fixed=TRUE)==TRUE)), pch=21, fill="#0072B2", size=3)+
    geom_tippoint(aes(
      subset=(grepl( "Frischella",label,fixed=TRUE)==TRUE)), pch=21, fill="#CC79A7", size=3)+
    geom_tippoint(aes(
      subset=(grepl("Bartonella",label,fixed=TRUE)==TRUE)), pch=21, fill="#D55E00", size=3) 
  
  


  # "#332288" "#6699CC"
  # [3] "#88CCEE" "#44AA99"
  # [5] "#117733" "#999933"
  # [7] "#DDCC77" "#661100"
  # [9] "#CC6677" "#AA4466"
  # [11] "#882255" "#AA4499"
  
  p2
  
  #browser()
}

## obligate symbionts
these_obligates <- c("Gilliamella",
                     "Lactobacillus",
                     "Snodgrassella",
                     "Commensalibacter",
                     "Bifidobacterium",
                     "Frischella",
                     "Bartonella")

# core.sp <- c("Rosenbergiella", "Pseudomonas", "Gilliamella",
#              "Lactobacillus", "Caulobacter", "Snodgrassella",
#              "Acinetobacter", "Corynebacterium", "Sphingomonas",
#              "Commensalibacter", "Methylobacterium",
#              "Massilia","Stenotrophomonas", "Bifidobacterium", "Frischella", "Bartonella")

### what i need to do:
# figure out what order to do the filtering steps
# featurs_site_metadata needs to be able to keep track of the tips that get renamed and if they are kept in the tree, attache sitecount metadata to the correct spot
# not sure how it is handling branches that are collapsed, they should be additive so that if fam1 was at 2 sites and fam2 was also at two sites, if they are the same sites the total number should be 2 but if they are different should be 4
# i think there are a few ways to do this
# what i will try next will be adding an index to keep track of the branch location and then after that dropping down the taxonomy and using the tip index to assign site count but surely this will be more complicated than that... :/ 

#tree.type <- 'rbcl'
meta_cols <- c('UniqueID', 'Family', 'Genus', 'Species', 'GenusSpecies', 'Sex', 'Site', 'Meadow')

meta <- spec.net %>%
  filter(Apidae == 1) %>%
    select(all_of(meta_cols), Apidae, starts_with('16s')) %>%
    na.omit() %>%
    select(!starts_with('16s')) 

comm_presabs <- as.data.frame(indiv.comm.16sR0) #load in the pres/abs table
comm_presabs[comm_presabs > 0] <- 1 #change all rel abund to 1
comm_presabs <- tibble::rownames_to_column(comm_presabs, "UniqueID") #make rownames (UniqueID) into column






# #apis
# # apis_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Apis", comm_presabs, apis_sites, all_levels=TRUE, do_collapse=FALSE, tree.type='16s')
# # apis_tree
# # 
# apis_tree2 <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Apis mellifera", genus.or.spp='Species', comm_presabs, apis_sites, all_levels=TRUE, do_collapse = TRUE)
# apis_tree2
# # 
# apis_tree_drop_s <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Apis mellifera", genus.or.spp='Species', comm_presabs, apis_sites, all_levels=FALSE, levels_to_drop=' s__', do_collapse=TRUE)
# apis_tree_drop_s
# # 
# apis_tree_drop_g <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Apis mellifera", genus.or.spp='Species', comm_presabs, apis_sites, all_levels=FALSE, levels_to_drop=c(' s__', ' g__'), do_collapse = TRUE)
# apis_tree_drop_g
# 
# 
# 
# #bombus centralis
# centralis_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus centralis", genus.or.spp='Species', comm_presabs, bombus_sites, do_collapse = TRUE)
# centralis_tree
# 
# centralis_tree_drop_s <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus centralis", genus.or.spp='Species', comm_presabs, bombus_sites, all_levels=FALSE, levels_to_drop=' s__', do_collapse = TRUE)
# centralis_tree_drop_s
# 
# centralis_tree_drop_g <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus centralis", genus.or.spp='Species', comm_presabs, bombus_sites, all_levels=FALSE, levels_to_drop=c(' s__', ' g__'), do_collapse = TRUE)
# centralis_tree_drop_g
# 
# #bombus huntii
# huntii_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus huntii", genus.or.spp='Species', comm_presabs, bombus_sites, do_collapse = TRUE)
# huntii_tree
# 
# huntii_tree_drop_s <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus huntii", genus.or.spp='Species', comm_presabs,  bombus_sites, all_levels=FALSE, levels_to_drop=' s__', do_collapse = TRUE)
# huntii_tree_drop_s
# 
# huntii_tree_drop_g <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus huntii", genus.or.spp='Species', comm_presabs, bombus_sites, all_levels=FALSE, levels_to_drop=c(' s__', ' g__'), do_collapse = TRUE)
# huntii_tree_drop_g
# 
# #bombus bifarius -- only found at one site!
# # bifarius_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus bifarius", comm_presabs, bombus_sites, do_collapse = TRUE)
# # bifarius_tree
# #
# # bifarius_tree_drop_s <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus bifarius", comm_presabs, bombus_sites, all_levels=FALSE, levels_to_drop=' s__', do_collapse = TRUE)
# # bifarius_tree_drop_s
# #
# # bifarius_tree_drop_g <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus bifarius", comm_presabs, bombus_sites, all_levels=FALSE, levels_to_drop=c(' s__', ' g__'), do_collapse = TRUE)
# # bifarius_tree_drop_g
# 
# 
# # melissodes confusus
# melissodes_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Melissodes confusus", genus.or.spp='Species', comm_presabs, bombus_sites, do_collapse = TRUE)
# melissodes_tree
# 
# melissodes_tree_drop_s <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Melissodes confusus", genus.or.spp='Species', comm_presabs, bombus_sites, all_levels=FALSE, levels_to_drop=' s__', do_collapse = TRUE)
# melissodes_tree_drop_s
# 
# melissodes_tree_drop_g <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Melissodes confusus", genus.or.spp='Species', comm_presabs, bombus_sites, all_levels=FALSE, levels_to_drop=c(' s__', ' g__'), do_collapse = TRUE)
# melissodes_tree_drop_g


## Genus trees
apis_tree2 <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Apis", genus.or.spp='Genus', comm_presabs, apis_sites, all_levels=TRUE, do_collapse = TRUE)
panelB <- apis_tree2 + labs(tag="B.")
panelB
# 
apis_tree_drop_s <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Apis", genus.or.spp='Genus', comm_presabs, apis_sites, all_levels=FALSE, levels_to_drop=' s__', do_collapse=TRUE)
apis_tree_drop_s
# 


# melissodes 
melissodes_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Melissodes", genus.or.spp='Genus', comm_presabs, bombus_sites, do_collapse = TRUE)
panelC <- melissodes_tree + labs(tag="C.")

melissodes_tree_drop_s <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Melissodes", genus.or.spp='Genus', comm_presabs, bombus_sites, all_levels=FALSE, levels_to_drop=' s__', do_collapse = TRUE)
melissodes_tree_drop_s


#bombus tree
bombus_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus", genus.or.spp='Genus', comm_presabs, bombus_sites, do_collapse = TRUE)
panelA <- bombus_tree + labs(tag="A.")

bombus_tree_drop_s <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus", genus.or.spp='Genus', comm_presabs,  bombus_sites, all_levels=FALSE, levels_to_drop=' s__', do_collapse = TRUE)
bombus_tree_drop_s



## make legend


# Load Packages 
library("ggplot2") 
library("grid") 
library("gridExtra") 
library("cowplot") 

core.sp <- c("Rosenbergiella", "Pseudomonas", "Gilliamella",
             "Lactobacillus", "Caulobacter", "Snodgrassella",
             "Acinetobacter", "Corynebacterium", "Sphingomonas",
             "Commensalibacter", "Methylobacterium",
             "Massilia","Stenotrophomonas", "Bifidobacterium", "Frischella", "Bartonella")

geom_tippoint(aes(
  subset=(grepl("Gilliamella",label,fixed=TRUE)==TRUE)), pch=21, fill="#E69F00", size=3)+
  geom_tippoint(aes(
    subset=(grepl("Lactobacillus",label,fixed=TRUE)==TRUE)), pch=21, fill="#56B4E9", size=3)+
  geom_tippoint(aes(
    subset=(grepl("Snodgrassella",label,fixed=TRUE)==TRUE)), pch=21, fill="#009E73", size=3)+
  geom_tippoint(aes(
    subset=(grepl("Commensalibacter",label,fixed=TRUE)==TRUE)), pch=21, fill="#F0E442", size=3)+
  geom_tippoint(aes(
    subset=(grepl( "Bifidobacterium",label,fixed=TRUE)==TRUE)), pch=21, fill="#0072B2", size=3)+
  geom_tippoint(aes(
    subset=(grepl( "Frischella",label,fixed=TRUE)==TRUE)), pch=21, fill="#CC79A7", size=3)+
  geom_tippoint(aes(
    subset=(grepl("Bartonella",label,fixed=TRUE)==TRUE)), pch=21, fill="#D55E00", size=3)

leg_col_order <- c("#D55E00", #Bartonella
                   "#0072B2", #Bifidobacterium
                   "#F0E442", #commensalibacter
                   "#CC79A7", #Frischella
                   "#E69F00", #gilliamella
                   "#56B4E9", #lactobacillus
                   "#009E73" )#snodgrassella



# Create a DataFrame 
data <- data.frame( 
  Xdata = rnorm(7),
  Ydata = rnorm(7), 
  Family = these_obligates,
  leg_color = leg_col_order)
                
    # "#E69F00", #gilliamella 5
    #             "#56B4E9", #Lactobacillus 6
    #             "#009E73", #Snodgrassella 7
    #             "#F0E442", #Commensalibacter3
    #             "#0072B2", #Bifidobacterium2
    #             "#CC79A7", #Frischella 4
    #             "#D55E00")) #Brtonella 1
    
    
    

# Create a Scatter Plot 
gplot <- ggplot(data, aes(Xdata, Ydata, color = Family)) +    
  geom_point(size = 7) +
  scale_color_manual(values=data$leg_color) +
  theme(legend.position='bottom') +
  labs(color='Bacteria Genus') +
  guides(colour = guide_legend(nrow = 1)) +  theme(legend.key=element_blank())

# Draw Only Legend without plot 
# Grab legend from gplot 
panelD  <- get_legend(gplot)                     
plot(get_legend(gplot) )

#make panel figure
grid.arrange(panelA,
             panelB,
             panelC,
             #panelD,
             ncol=3)
             #layout_matrix = cbind(c(1,4), c(2,4), c(3,4)))

