
## tutorial: https://bioconductor.org/packages/devel/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html



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
library(viridis)

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
  filter(Apidae == 1) #%>%
  # mutate(Site = factor(Site, levels=c("JC", ## ordered by latitude north to south
  #                                     "SM",
  #                                     "SC",
  #                                     "MM",
  #                                     "HM",
  #                                     "PL",
  #                                     "CH")))


##can probs figure out a way to do this automatically but was getting frustrated so hard coded this
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

#####

##import community presence/absence file
setwd("../../skyIslands_saved")

comm_presabs <- as.data.frame(indiv.comm.16sR0) #load in the pres/abs table
comm_presabs[comm_presabs > 0] <- 1 #change all rel abund to 1
comm_presabs <- tibble::rownames_to_column(comm_presabs, "UniqueID") #make rownames (UniqueID) into column

#######################################


##function that filters a tree object to a certain sample bee genus, generates a presence absence heatmap for whether that tree tip was found 
## at any given site, then plots the tree with the appended heatmap


##probs want to make sure each site in each individual graph is colored uniformly for all graphs 

phylotree_heatmap_byGenus <- function(tree.object, metadata, this.genus, presAbsTable, site.order){
  genus_ids <- metadata %>%
    filter(Genus==this.genus) %>%
    select(UniqueID)
  
  my_sites <- unique(metadata$Site[metadata$Genus==this.genus])

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

  
  matched_presabs <- match_shared_tiplabels(gentree, presAbsTable)
  
  matched_pres_meta <- match_shared_ID(matched_presabs, metadata)
  
  matched_id <- matched_pres_meta$UniqueID
  row.names(matched_pres_meta) <- matched_id
  
  meta_match_sites <- match_shared_ID(metadata, matched_pres_meta) %>%
    select(UniqueID, Site, Genus) %>%
    mutate(Site = factor(Site)) %>%
    filter(Genus==this.genus) %>%
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
    filter(Site_present > 0) %>%
    mutate(Site = factor(Site, levels=site.order))
  
  features_site_metadata <- match_shared_ID(matched_pres_meta, meta_match_sites) %>%
    right_join(meta_match_sites, by='UniqueID') %>%
    pivot_longer(cols = starts_with('16s'), names_to = 'bacteria', values_to = 'bact_pres') %>%
    group_by(bacteria) %>%
    filter(bact_pres == 1) %>%
    select(!bact_pres) %>%
    relocate(bacteria)
  
  
  p <- ggtree(gentree, layout='rectangular') 
  p
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
        text.size=2,
        vjust=-110,
        #text.angle=-45
      ),
      show.legend=FALSE) +
    scale_fill_viridis(option="plasma", discrete=TRUE) +
    ggtitle(this.genus)
  p2
}

apis_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Apis", comm_presabs, apis_sites)
apis_tree


bombus_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Bombus", comm_presabs, bombus_sites)
bombus_tree

anthophora_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Anthophora",comm_presabs, anthophora_sites)
anthophora_tree

megachile_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Megachile", comm_presabs, megachile_sites)
megachile_tree

##########################

##getting clade labels ready

##used 

apis_table <- apis_tree%>% as.treedata %>% as_tibble


apis_with_clades <- apis_tree + 
  geom_cladelab(node=209, label="Bifidobacteriaceae", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .1, textcolor='red', barcolor='red') +
  geom_cladelab(node=236, label="Lactobacillaceae", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .1, textcolor='orange', barcolor='orange') +
  geom_cladelab(node=277, label="Neisseriaceae", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .1, textcolor='gold', barcolor='gold') +
  geom_cladelab(node=285, label="Orbaceae", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .1, textcolor='green', barcolor='green')  +
  geom_cladelab(node=312, label="Enterobacteriaceae", angle=270, offset=.6, hjust='center', align=TRUE, offset.text = .1, textcolor='blue', barcolor='blue')+
  geom_cladelab(node=357, label="Bartonella", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .2, textcolor='violet', barcolor='violet') +
  geom_cladelab(node=364, label="Acetobacteriaceae", angle=270, offset=.6, hjust='center', align=TRUE, offset.text = .3, textcolor='purple', barcolor='purple') 
  
apis_with_clades 

###
bombus_table <- bombus_tree%>% as.treedata %>% as_tibble


bombus_with_clades <- bombus_tree + 
  geom_cladelab(node=679, label="Bifidobacteriaceae", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .1, textcolor='red', barcolor='red') +
  #geom_cladelab(node=236, label="Lactobacillaceae", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .1, textcolor='orange', barcolor='orange') +
  #geom_cladelab(node=277, label="Neisseriaceae", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .1, textcolor='gold', barcolor='gold') +
  #geom_cladelab(node=285, label="Orbaceae", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .1, textcolor='green', barcolor='green')  +
  #geom_cladelab(node=312, label="Enterobacteriaceae", angle=270, offset=.6, hjust='center', align=TRUE, offset.text = .1, textcolor='blue', barcolor='blue')+
  geom_cladelab(node=1167, label="Bartonella", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .2, textcolor='violet', barcolor='violet') +
  #geom_cladelab(node=364, label="Acetobacteriaceae", angle=270, offset=.6, hjust='center', align=TRUE, offset.text = .3, textcolor='purple', barcolor='purple') 

bombus_with_clades

##maybe use geom_strip for non monophyletic groups --
##can i write a function that searches tip labels for D_4 or D_5 string pattern for correct family or genus of interest,
##returns the node numbers or tip numbers/labels based on a list of correct labels
## probs need a for loop to do this, input will be a list of labels
## 1. for loop
## 2. 1:length of list of clades
## 3. first need to search tip labels for pattern matches
## 4. return either tipnumbers or node numbers 
## then need to add geom)cladelab with the node
## if not monophyletic then need to add a second bar based on the tip numbers/labels with matched color

####################### copying down apis function as practice example



genus_ids <- meta %>%
  filter(Genus=='Apis') %>%
  select(UniqueID)

my_sites <- unique(meta$Site[meta$Genus=='Apis'])

trimmed_tree <- prune_samples(sample_names(physeq16sR0) %in% genus_ids$UniqueID, physeq16sR0)

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


true_tips <- grepl('D_4__Orbaceae', gentree$tip.label) #boolean to determine which tip labels match the fam of interest

fam_tips <- gentree$tip.label[true_tips] #filter to just those labels

tree_table <- p %>% as.treedata %>% as_tibble

family_nodes <- tree_table$node[tree_table$label %in% fam_tips == TRUE]

p1 <- p + 
  geom_tippoint(aes(subset=p$data$node %in% family_nodes), color='red', size=0.5)
p1

##########
## 1-28 struggling to functionize this... maybe just add one by one in 
## the function? I think i would like the tips to be colored 
## and there to be bars showing labs... maybe bring back code
## from the other day that kinda clunkily did this?

## pasted below

##really what i want is:
## add to the function the ability to
## input a list of bacteria tips labels of interest
## for each, color only tips with labels that include that family
## need to add a legend or clade label
## clade label bar can be imperfect since we are 
## coloring actual tips -- hope this will be a workaround
## for weirdness with nonmonophyletic taxa

## have found the correct nodes need to figure out how to segment them
## want to:
# 1. make a list of every sequential combo (index1-index2, index2-index3 etc)

#Create nested list (3 lists inside)
my_nested_list <- list(taxa1=apis_nodes[1:(length(apis_nodes)-1)],
                       taxa2=apis_nodes[2:length(apis_nodes)])

# Convert nested list to the dataframe by columns
node_comparisons <- as.data.frame(do.call(cbind, my_nested_list))
node_comparisons

#now have appropriate list of node segments, need to 
#find a way to draw segments between each or otherwise figure
#out how to annotate the polyphyletic groups 

##this below for loop isnt working properly yet

# for (i in seq_along(node_comparisons)){
#   apis_with_clades <- apis_tree + 
#     geom_strip(taxa1=node_comparisons$taxa1[i],
#                taxa2=node_comparisons$taxa2[i],
#                label="")
# }
# apis_with_clades

# apis_with_clades <- apis_tree + 
#   geom_strip(apis_nodes, label="Lactobacillaceae", angle=270, hjust='center', offset=.6, align=TRUE, offset.text = .1, textcolor='red', barcolor='red')
