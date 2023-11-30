
## tutorial: https://bioconductor.org/packages/devel/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html
rm(list=ls())

#####
setwd('~/Dropbox (University of Oregon)/skyIslands/')


source('dataPrep/src/trees_init.R')
load('presAbsTable.Rdata')
load('spec_RBCL_16s.Rdata')
load('physeq16s.Rdata')
##import community presence/absence file
setwd("../../skyIslands_saved")

library(pals)


#######################################


##function that filters a tree object to a certain sample bee genus, generates a presence absence heatmap for whether that tree tip was found 
## at any given site, then plots the tree with the appended heatmap


##probs want to make sure each site in each individual graph is colored uniformly for all graphs 

phylotree_heatmap_byGenus <- function(tree.object, metadata, this.species, presAbsTable, site.order, all_levels=TRUE, levels_to_drop, clade_names=NULL, do_collapse=FALSE){
  #filter to include just the unique IDs in the specified genus
  sp_ids <- metadata %>%
    filter(GenusSpecies==this.species) %>%
    select(UniqueID)
  
  #pull out all sites that include the specified genus
  my_sites <- unique(metadata$Site[metadata$GenusSpecies==this.species])

  #remove tips from the tree that are not in the list of unique IDs in the specified genus
  trimmed_tree <- prune_samples(sample_names(tree.object) %in% sp_ids$UniqueID, tree.object)
  
  #remove taxa from the tree where there were less than zero observations
  pruned_tree <- prune_taxa(taxa_sums(trimmed_tree) > 0, trimmed_tree)
  
  #read in the taxonomic info for each feature
  feature.2.tax.16s <-
    read.table("SI_pipeline/R2018/2023_sequence_results_raw/merged/16s/taxonomy.tsv", sep="\t",
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
      rest_of_label <- ' to genus'
    } else if(final_level == ' g__'){
      rest_of_label <- ' to family'
    }
    this_level <- paste(": Collapsed", rest_of_label)
  } else {this_level <- ': Full Tree'}

  
  matched_presabs <- match_shared_tiplabels(gentree, presAbsTable)
  
  matched_pres_meta <- match_shared_ID(matched_presabs, metadata)
  
  matched_id <- matched_pres_meta$UniqueID
  row.names(matched_pres_meta) <- matched_id
  
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
    mutate(SiteCount = as.numeric(n_distinct(Site))) 
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
    scale_fill_gradient(high = "black", low ="lightgrey") +
    ggtitle(paste(this.species, this_level)) +
    new_scale_fill() + 
    geom_tippoint(aes(
      subset=(grepl('Acetobacteraceae',label,fixed=TRUE)==TRUE)), color="#604E97", size=3) +
    geom_tippoint(aes(
      subset=(grepl('Bifidobacteriaceae',label,fixed=TRUE)==TRUE)), color="#F3C300", size=3)+
    geom_tippoint(aes(
      subset=(grepl('Enterobacteriaceae',label,fixed=TRUE)==TRUE)), color="#B3446C", size=3)+
    geom_tippoint(aes(
      subset=(grepl('Erwiniaceae',label,fixed=TRUE)==TRUE)), color="#F38400", size=3)+
    geom_tippoint(aes(
      subset=(grepl('Hafniaceae',label,fixed=TRUE)==TRUE)), color="#A1CAF1", size=3)+
    geom_tippoint(aes(
      subset=(grepl('Lactobacillaceae',label,fixed=TRUE)==TRUE)), color="#BE0032", size=3)+
    geom_tippoint(aes(
      subset=(grepl('Leuconostocaceae',label,fixed=TRUE)==TRUE)), color="#C2B280", size=3)+
    geom_tippoint(aes(
      subset=(grepl('Moraxellaceae',label,fixed=TRUE)==TRUE)), color="#848482", size=3)+
    geom_tippoint(aes(
      subset=(grepl('Neisseriaceae',label,fixed=TRUE)==TRUE)), color="#008856", size=3)+
    geom_tippoint(aes(
      subset=(grepl('Orbaceae',label,fixed=TRUE)==TRUE)), color="#E68FAC", size=3)+
    geom_tippoint(aes(
      subset=(grepl('Streptococcaceae',label,fixed=TRUE)==TRUE)), color="#0067A5", size=3)+
    geom_tippoint(aes(
      subset=(grepl('Yersiniaceae',label,fixed=TRUE)==TRUE)), color="#F99379", size=3)
  


  # "#332288" "#6699CC"
  # [3] "#88CCEE" "#44AA99"
  # [5] "#117733" "#999933"
  # [7] "#DDCC77" "#661100"
  # [9] "#CC6677" "#AA4466"
  # [11] "#882255" "#AA4499"
  
  p2
  
  #browser()
}

### what i need to do:
# figure out what order to do the filtering steps
# featurs_site_metadata needs to be able to keep track of the tips that get renamed and if they are kept in the tree, attache sitecount metadata to the correct spot
# not sure how it is handling branches that are collapsed, they should be additive so that if fam1 was at 2 sites and fam2 was also at two sites, if they are the same sites the total number should be 2 but if they are different should be 4
# i think there are a few ways to do this
# what i will try next will be adding an index to keep track of the branch location and then after that dropping down the taxonomy and using the tip index to assign site count but surely this will be more complicated than that... :/ 

#tree.type <- 'rbcl'

meta <- spec.net %>%
  filter(Apidae == 1) %>%
    select(all_of(meta_cols), Apidae, starts_with('16s')) %>%
    na.omit() %>%
    select(!starts_with('16s')) 

comm_presabs <- as.data.frame(indiv.comm.16sR0) #load in the pres/abs table
comm_presabs[comm_presabs > 0] <- 1 #change all rel abund to 1
comm_presabs <- tibble::rownames_to_column(comm_presabs, "UniqueID") #make rownames (UniqueID) into column






#apis
# apis_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Apis", comm_presabs, apis_sites, all_levels=TRUE, do_collapse=FALSE, tree.type='16s')
# apis_tree
# 
apis_tree2 <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Apis mellifera", comm_presabs, apis_sites, all_levels=TRUE, do_collapse = TRUE)
apis_tree2
# 
apis_tree_drop_s <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Apis mellifera", comm_presabs, apis_sites, all_levels=FALSE, levels_to_drop=' s__', do_collapse=TRUE)
apis_tree_drop_s
# 
apis_tree_drop_g <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Apis mellifera", comm_presabs, apis_sites, all_levels=FALSE, levels_to_drop=c(' s__', ' g__'), do_collapse = TRUE)
apis_tree_drop_g



#bombus centralis
centralis_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus centralis", comm_presabs, bombus_sites, do_collapse = TRUE)
centralis_tree

centralis_tree_drop_s <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus centralis", comm_presabs, bombus_sites, all_levels=FALSE, levels_to_drop=' s__', do_collapse = TRUE)
centralis_tree_drop_s

centralis_tree_drop_g <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus centralis", comm_presabs, bombus_sites, all_levels=FALSE, levels_to_drop=c(' s__', ' g__'), do_collapse = TRUE)
centralis_tree_drop_g

#bombus huntii
huntii_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus huntii", comm_presabs, bombus_sites, do_collapse = TRUE)
huntii_tree

huntii_tree_drop_s <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus huntii", comm_presabs, bombus_sites, all_levels=FALSE, levels_to_drop=' s__', do_collapse = TRUE)
huntii_tree_drop_s

huntii_tree_drop_g <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus huntii", comm_presabs, bombus_sites, all_levels=FALSE, levels_to_drop=c(' s__', ' g__'), do_collapse = TRUE)
huntii_tree_drop_g

#bombus bifarius -- only found at one site!
# bifarius_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus bifarius", comm_presabs, bombus_sites, do_collapse = TRUE)
# bifarius_tree
#
# bifarius_tree_drop_s <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus bifarius", comm_presabs, bombus_sites, all_levels=FALSE, levels_to_drop=' s__', do_collapse = TRUE)
# bifarius_tree_drop_s
#
# bifarius_tree_drop_g <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus bifarius", comm_presabs, bombus_sites, all_levels=FALSE, levels_to_drop=c(' s__', ' g__'), do_collapse = TRUE)
# bifarius_tree_drop_g


# melissodes confusus
melissodes_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Melissodes confusus", comm_presabs, bombus_sites, do_collapse = TRUE)
melissodes_tree

melissodes_tree_drop_s <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Melissodes confusus", comm_presabs, bombus_sites, all_levels=FALSE, levels_to_drop=' s__', do_collapse = TRUE)
melissodes_tree_drop_s

melissodes_tree_drop_g <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Melissodes confusus", comm_presabs, bombus_sites, all_levels=FALSE, levels_to_drop=c(' s__', ' g__'), do_collapse = TRUE)
melissodes_tree_drop_g

## make legend


# Load Packages 
library("ggplot2") 
library("grid") 
library("gridExtra") 
library("cowplot") 

geom_tippoint(aes(
  subset=(grepl('Acetobacteraceae',label,fixed=TRUE)==TRUE)), color="#604E97", size=3) +
  geom_tippoint(aes(
    subset=(grepl('Bifidobacteriaceae',label,fixed=TRUE)==TRUE)), color="#F3C300", size=3)+
  geom_tippoint(aes(
    subset=(grepl('Enterobacteriaceae',label,fixed=TRUE)==TRUE)), color="#B3446C", size=3)+
  geom_tippoint(aes(
    subset=(grepl('Erwiniaceae',label,fixed=TRUE)==TRUE)), color="#F38400", size=3)+
  geom_tippoint(aes(
    subset=(grepl('Hafniaceae',label,fixed=TRUE)==TRUE)), color="#A1CAF1", size=3)+
  geom_tippoint(aes(
    subset=(grepl('Lactobacillaceae',label,fixed=TRUE)==TRUE)), color="#BE0032", size=3)+
  geom_tippoint(aes(
    subset=(grepl('Leuconostocaceae',label,fixed=TRUE)==TRUE)), color="#C2B280", size=3)+
  geom_tippoint(aes(
    subset=(grepl('Moraxellaceae',label,fixed=TRUE)==TRUE)), color="#848482", size=3)+
  geom_tippoint(aes(
    subset=(grepl('Neisseriaceae',label,fixed=TRUE)==TRUE)), color="#008856", size=3)+
  geom_tippoint(aes(
    subset=(grepl('Orbaceae',label,fixed=TRUE)==TRUE)), color="#E68FAC", size=3)+
  geom_tippoint(aes(
    subset=(grepl('Streptococcaceae',label,fixed=TRUE)==TRUE)), color="#0067A5", size=3)+
  geom_tippoint(aes(
    subset=(grepl('Yersiniaceae',label,fixed=TRUE)==TRUE)), color="#F99379", size=3)




# Create a DataFrame 
data <- data.frame( 
  Xdata = rnorm(12),
  Ydata = rnorm(12), 
  Family = c('Acetobacteraceae',
                 'Bifidobacteriaceae',
                 'Enterobacteriaceae',
                 'Erwiniaceae',
                 'Hafniaceae', 
                 'Lactobacillaceae',
                 'Leuconostocaceae',
                 'Moraxellaceae',  
                 'Neisseriaceae',
                 'Orbaceae',
                 'Streptococcaceae',
                 'Yersiniaceae'
                 ),
  leg_color = c("#604E97",
                "#F3C300",
                "#B3446C",
                "#F38400",
                "#A1CAF1",
                "#BE0032",
                "#C2B280",
                "#848482",
                "#008856",
                "#E68FAC",
                "#0067A5",
                "#F99379"
                )) 

# Create a Scatter Plot 
gplot <- ggplot(data, aes(Xdata, Ydata, color = Family)) +    
  geom_point(size = 7) + scale_color_manual(values=data$leg_color)

# Draw Only Legend without plot 
# Grab legend from gplot 
legend <- get_legend(gplot)                     

# Create new plot window 
grid.newpage()                               

# Draw Only legend  
grid.draw(legend) 
# ##########################
# 
# ##getting clade labels ready
# 
# ##used 
# 
# apis_table <- apis_tree%>% as.treedata %>% as_tibble
# 
# 
# apis_with_clades <- apis_tree + 
#   geom_cladelab(node=209, label="Bifidobacteriaceae", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .1, textcolor='red', barcolor='red') +
#   geom_cladelab(node=236, label="Lactobacillaceae", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .1, textcolor='orange', barcolor='orange') +
#   geom_cladelab(node=277, label="Neisseriaceae", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .1, textcolor='gold', barcolor='gold') +
#   geom_cladelab(node=285, label="Orbaceae", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .1, textcolor='green', barcolor='green')  +
#   geom_cladelab(node=312, label="Enterobacteriaceae", angle=270, offset=.6, hjust='center', align=TRUE, offset.text = .1, textcolor='blue', barcolor='blue')+
#   geom_cladelab(node=357, label="Bartonella", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .2, textcolor='violet', barcolor='violet') +
#   geom_cladelab(node=364, label="Acetobacteriaceae", angle=270, offset=.6, hjust='center', align=TRUE, offset.text = .3, textcolor='purple', barcolor='purple') 
#   
# apis_with_clades 
# 
# ###
# bombus_table <- bombus_tree%>% as.treedata %>% as_tibble
# 
# 
# bombus_with_clades <- bombus_tree + 
#   geom_cladelab(node=679, label="Bifidobacteriaceae", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .1, textcolor='red', barcolor='red') +
#   #geom_cladelab(node=236, label="Lactobacillaceae", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .1, textcolor='orange', barcolor='orange') +
#   #geom_cladelab(node=277, label="Neisseriaceae", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .1, textcolor='gold', barcolor='gold') +
#   #geom_cladelab(node=285, label="Orbaceae", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .1, textcolor='green', barcolor='green')  +
#   #geom_cladelab(node=312, label="Enterobacteriaceae", angle=270, offset=.6, hjust='center', align=TRUE, offset.text = .1, textcolor='blue', barcolor='blue')+
#   geom_cladelab(node=1167, label="Bartonella", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .2, textcolor='violet', barcolor='violet') +
#   #geom_cladelab(node=364, label="Acetobacteriaceae", angle=270, offset=.6, hjust='center', align=TRUE, offset.text = .3, textcolor='purple', barcolor='purple') 
# 
# bombus_with_clades
# 
# ##maybe use geom_strip for non monophyletic groups --
# ##can i write a function that searches tip labels for D_4 or D_5 string pattern for correct family or genus of interest,
# ##returns the node numbers or tip numbers/labels based on a list of correct labels
# ## probs need a for loop to do this, input will be a list of labels
# ## 1. for loop
# ## 2. 1:length of list of clades
# ## 3. first need to search tip labels for pattern matches
# ## 4. return either tipnumbers or node numbers 
# ## then need to add geom)cladelab with the node
# ## if not monophyletic then need to add a second bar based on the tip numbers/labels with matched color
# 
# ####################### copying down apis function as practice example
# 
# 
# genus_ids <- meta %>%
#   filter(Genus=='Apis') %>%
#   select(UniqueID)
# 
# my_sites <- unique(meta$Site[meta$Genus=='Apis'])
# 
# trimmed_tree <- prune_samples(sample_names(physeq16sR0) %in% genus_ids$UniqueID, physeq16sR0)
# 
# pruned_tree <- prune_taxa(taxa_sums(trimmed_tree) > 0, trimmed_tree)
# 
# 
# feature.2.tax.16s <-
#   read.table("SI_pipeline/merged/16s/taxonomy16s.txt", sep="\t",
#              header=TRUE)
# 
# feature.2.tax.16s$Taxon <- paste("16s", feature.2.tax.16s$Taxon, sep=':')
# 
# # ## convert to a phylo class which is more useful downstream
# gentree <- phy_tree(pruned_tree, errorIfNULL=TRUE)
# 
# ## match the tip labs to the table with feature ID and Taxon
# gentree$tip.label  <-  feature.2.tax.16s$Taxon[match(gentree$tip.label,
#                                                      feature.2.tax.16s$Feature.ID)]
# 
# p <- ggtree(gentree, layout='rectangular') 
# p
# 
# 
# ##########
# ## 1-28 struggling to functionize this... maybe just add one by one in 
# ## the function? I think i would like the tips to be colored 
# ## and there to be bars showing labs... maybe bring back code
# ## from the other day that kinda clunkily did this?
# 
# ## pasted below
# 
# ## 1-30 works below by hand... struggling to turn into function, probs will have to do by
# ## hand in the interest of time
# 
# node_list <- list() #initialize list for nodes
# 
# for (fam.name in bact_fam_list){
#   
#   true_tips <- grepl(fam.name, gentree$tip.label) #boolean to determine which tip labels match the fam of interest
#   
#   fam_tips <- gentree$tip.label[true_tips] #filter to just those labels
#   
#   tree_table <- p %>% as.treedata %>% as_tibble
#   
#   family_nodes <- tree_table$node[tree_table$label %in% fam_tips == TRUE]
#   
#   family_nodes <- list(family_nodes)
#   
#   node_list <- append(node_list, family_nodes)
#   
# }
# 
# node_list
# 
# for (i in 1:length(node_list)){
#   
#   these_nodes <- p$data$node %in% unlist(node_list[[i]])
#   
#   tip_plot <- p + geom_tippoint(aes(subset=these_nodes),
#                                 color=randomColor(),
#                                 size=0.7)
#   
#   tip_plot
#   
#   p <- tip_plot
# }
# 
# 
# ##really what i want is:
# ## add to the function the ability to
# ## input a list of bacteria tips labels of interest
# ## for each, color only tips with labels that include that family
# ## need to add a legend or clade label
# ## clade label bar can be imperfect since we are 
# ## coloring actual tips -- hope this will be a workaround
# ## for weirdness with nonmonophyletic taxa
# 
# ## have found the correct nodes need to figure out how to segment them
# ## want to:
# # 1. make a list of every sequential combo (index1-index2, index2-index3 etc)
# 
# #Create nested list (3 lists inside)
# my_nested_list <- list(taxa1=apis_nodes[1:(length(apis_nodes)-1)],
#                        taxa2=apis_nodes[2:length(apis_nodes)])
# 
# # Convert nested list to the dataframe by columns
# node_comparisons <- as.data.frame(do.call(cbind, my_nested_list))
# node_comparisons
# 
# #now have appropriate list of node segments, need to 
# #find a way to draw segments between each or otherwise figure
# #out how to annotate the polyphyletic groups 
# 
# ##this below for loop isnt working properly yet
# 
# # for (i in seq_along(node_comparisons)){
# #   apis_with_clades <- apis_tree + 
# #     geom_strip(taxa1=node_comparisons$taxa1[i],
# #                taxa2=node_comparisons$taxa2[i],
# #                label="")
# # }
# # apis_with_clades
# 
# # apis_with_clades <- apis_tree + 
# #   geom_strip(apis_nodes, label="Lactobacillaceae", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .1, textcolor='red', barcolor='red')
# 
# 
# ##################################################
# 
# #phylotree_heatmap_byGenus(physeq16sR0, meta, "Apis", comm_presabs, apis_sites, all_levels=TRUE, clade_names='D_4__Lactobacillaceae')
# #function(physeq16sR0, metadata, this.genus, comm_presabs, apis_sites, all_levels=TRUE, levels_to_drop, clade_names=NULL){
# 
# 
# all_levels=TRUE
# 
# #figuring out adding triangle clade labs
# genus_ids <- meta %>%
#   filter(Genus=="Apis") %>%
#   select(UniqueID)
# 
# my_sites <- unique(meta$Site[meta$Genus=="Apis"])
# 
# trimmed_tree <- prune_samples(sample_names(physeq16sR0) %in% genus_ids$UniqueID, physeq16sR0)
# 
# pruned_tree <- prune_taxa(taxa_sums(trimmed_tree) > 0, trimmed_tree)
# 
# 
# feature.2.tax.16s <-
#   read.table("SI_pipeline/merged/16s/taxonomy16s.txt", sep="\t",
#              header=TRUE)
# 
# feature.2.tax.16s$Taxon <- paste("16s", feature.2.tax.16s$Taxon, sep=':')
# 
# # ## convert to a phylo class which is more useful downstream
# gentree <- phy_tree(pruned_tree, errorIfNULL=TRUE)
# 
# ## match the tip labs to the table with feature ID and Taxon
# gentree$tip.label  <-  feature.2.tax.16s$Taxon[match(gentree$tip.label,
#                                                      feature.2.tax.16s$Feature.ID)]
# 
# 
# # Identify tips with labels exactly matching '16s:D_0__Bacteria'
# matching_tips <- grep('^16s:D_0__Bacteria$', gentree$tip.label)
# 
# # Drop the matching tips
# gentree <- drop.tip(gentree, matching_tips)
# orig_len <- length(gentree$tip.label)
# print(paste('original # tips', orig_len))
# 
# 
# if(all_levels==FALSE){
#   for (level in levels_to_drop) {
#     tips_to_drop <- grep(level, gentree$tip.label)
#     print(paste('tip drop len', length(tips_to_drop)))
#     if (length(tips_to_drop) > 0) {
#       gentree <- drop.tip(gentree, tips_to_drop)
#       new_len <- length(gentree$tip.label)
#       percent_drop <- (length(tips_to_drop)/orig_len)*100
#       print(paste(level, "percent tips dropped", percent_drop))
#       final_level <- level
#     }
#   }
# }
# 
# if(all_levels==FALSE){
#   this_level <- paste(": Collapsed", final_level)
# } else {this_level <- ': Full Tree'}
# 
# 
# matched_presabs <- match_shared_tiplabels(gentree, comm_presabs)
# 
# matched_pres_meta <- match_shared_ID(matched_presabs, meta)
# 
# matched_id <- matched_pres_meta$UniqueID
# row.names(matched_pres_meta) <- matched_id
# 
# meta_match_sites <- match_shared_ID(meta, matched_pres_meta) %>%
#   select(UniqueID, Site, Genus) %>%
#   mutate(Site = factor(Site)) %>%
#   filter(Genus=="Apis") %>%
#   select(!Genus) %>%
#   group_by(UniqueID, Site) %>%
#   count() %>%
#   pivot_wider(names_from=Site,
#               values_from = n,
#               names_expand = TRUE,
#               id_expand=TRUE) %>%
#   pivot_longer(cols=2:length(colnames(.)),
#                names_to='Site',
#                values_to='Site_present') %>%
#   filter(Site_present > 0) #%>%
# #mutate(Site = factor(Site, levels=apis_sites))
# #browser()
# 
# features_site_metadata <- match_shared_ID(matched_pres_meta, meta_match_sites) %>%
#   right_join(meta_match_sites, by='UniqueID') %>%
#   pivot_longer(cols = starts_with('16s'), names_to = 'bacteria', values_to = 'bact_pres') %>%
#   group_by(bacteria) %>%
#   filter(bact_pres == 1) %>%
#   select(!bact_pres) %>%
#   relocate(bacteria) %>%
#   group_by(bacteria) %>%
#   add_count(Site, name="n_individuals") %>%
#   mutate(SiteCount = as.numeric(n_distinct(Site))) #%>%
# #select(!Site)
# 
# #browser()
# 
# p <- ggtree(gentree, layout='rectangular') 
# p
# 
# p2 <- p +
#   new_scale_fill() + 
#   geom_tiplab(align=TRUE, size=1.5) + 
#   coord_cartesian(clip="off") +
#   geom_fruit(
#     data=features_site_metadata,
#     geom=geom_tile,
#     pwidth=0.1,
#     offset=0.7,
#     mapping=aes(y=bacteria,
#                 #x=SiteCount,
#                 fill=SiteCount, width=0.1),
#     show.legend=TRUE) +
#   scale_fill_viridis(option="rocket", discrete=FALSE, direction=-1) +
#   ggtitle(paste("Apis", this_level)) 
# p2
# 
# 
# # Step 1: Find the best MRCA that matches to the keywords or search patten
# 
# groups <- unique(gentree$tip.label)
# 
# collapse_identical_tips <- function(phy,tip_label){
#   #matching_tips is initialized with the indices of tips in the phylogenetic tree (phy) whose labels match the provided tip_label. The function identifies all tips with the same label.
#   matching_tips <- which(phy$tip.label==tip_label)
#   nt <- length(phy$tip.label) # number of tips in tree
#   nm <- length(matching_tips) # Number of tips matching the label
#   keep <- numeric(nm) #keep is initialized as a numeric vector of length nm. It is used to keep track of which tips should be retained (1) and which tips should be dropped (0) in the new tree.
#   
#   #The while loop iterates through the indices of matching_tips to determine which tips to keep and which to drop.
#   cur_tip <- 1
#   #Inside the loop, the variable cur_tip is the current tip being considered, and next_tip is the tip immediately after cur_tip.
#   while(cur_tip<=nm){
#     if(cur_tip == nm){
#       keep[cur_tip] <- 1
#       break
#     }
#     next_tip <- cur_tip + 1
#     #mrca_ (most recent common ancestor) is calculated using the getMRCA function for the tips identified by matching_tips[cur_tip] and matching_tips[next_tip]. This helps find the common ancestor of the current tip and the next tip.
#     mrca_ <- getMRCA(phy,c(matching_tips[cur_tip],matching_tips[next_tip]))
#     #descendants contains the indices of all descendants of the common ancestor, which includes both tips and internal nodes.
#     descendants <- getDescendants(phy, mrca_)
#     #descendant_tips is calculated to include only those indices from descendants that correspond to actual tips (i.e., indices less than or equal to nt).
#     descendant_tips <- descendants[descendants<=nt]
#     #The function checks if all descendant_tips are in the list of matching_tips. If they are, it means all these tips can be collapsed into a single branch, and they are marked to be kept.
#     #The variable keep is updated accordingly, and cur_tip is advanced to skip the tips that have been collapsed.
#     if(all(descendant_tips %in% matching_tips)){
#       keep[cur_tip] <- 1
#       cur_tip <- cur_tip + length(descendant_tips)
#     }else{ #If not all descendant_tips are in the list of matching_tips, it means that not all tips can be collapsed, so the current tip is marked to be kept, and cur_tip is incremented by 1.
#       keep[cur_tip] <- 1
#       cur_tip <- cur_tip + 1
#     }
#   }
#   #After the loop completes, the to_drop variable contains the indices of the tips that need to be dropped to collapse identical labels.
#   to_drop <- matching_tips[!keep]
#   #The function then creates a new phylogenetic tree (new_phy) by using the drop.tip function to remove the tips identified in to_drop.
#   new_phy <- drop.tip(phy,to_drop)
#   #browser()
#   #Finally, the new phylogenetic tree is returned as the output of the function.
#   return(new_phy)
# }
# 
# for (this.group in groups){
#   gentree <- collapse_identical_tips(gentree, this.group)
# }
# 
