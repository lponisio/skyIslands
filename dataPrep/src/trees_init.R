

## packages
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ggtreeExtra")


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

## working dir

wdpath <- 'C:/Users/rah10/Dropbox (University of Oregon)/skyIslands/data'
setwd(wdpath)

#rm(list=ls())
## Data imports

spec16s <- read.csv('spec_RBCL_16s.csv') %>%
  filter(Apidae==1) 

names(spec16s) <- gsub(x = names(spec16s), pattern = "X1", replacement = "1")  
names(spec16s) <- gsub(x = names(spec16s), pattern = "RBCL.", replacement = "RBCL:")  

meta_cols <- c('UniqueID', 'Family', 'Genus', 'Species', 'Sex', 'Site', 'Meadow')

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
