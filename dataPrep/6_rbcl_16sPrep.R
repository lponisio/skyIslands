## from this
## tutorial:https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121

## tutorial also have some good pcoa plots and heatmaps, which arent
## on here
#dir.bombus <- '/Volumes/bombus/Dropbox (University of Oregon)'
rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path

setwd(local.path)
setwd('skyIslands_saved')



if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
##BiocManager::install("TreeSummarizedExperiment")
library(TreeSummarizedExperiment)
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)

library(tidyr)
library(dplyr)
library(bipartite)
library(phyloseq)
library(TreeTools)
library(devtools)
library(ape)
library(picante)

source("../skyIslands/dataPrep/src/misc.R")

# 2018 Samples!!
## reading artifacts
qza.16s.path  <- "SI_pipeline/merged/16s/final"
qza.rbcl.path  <- "SI_pipeline/merged/RBCL/final"

# 16s
weightedUF16sqzaR0 <- read_qza(file.path(qza.16s.path,
                      'core_metrics16s/weighted_unifrac_distance_matrix.qza'))

unweightedUF16sqzaR0 <- read_qza(file.path(qza.16s.path,
                      'core_metrics16s/unweighted_unifrac_distance_matrix.qza'))

# 16s
wphylo.dist.16sR0 <- weightedUF16sqzaR0$data
phylo.dist.16sR0 <-unweightedUF16sqzaR0$data

## Note, when taxonomy is imported, a single string is returned along
## with a confidence score.  For many analysis we will want to break
## up this string and for that purpose the parse_taxonomy() function
## is provided:

taxonomy16sR0 <- read_qza(
    file.path(qza.16s.path,
              "core_metrics16s/rarefied_table.qza"))

taxonomy16sR0 <- taxonomy16sR0$data

## ## 16s R0 phylogeny

physeq16sR0 <- qza_to_phyloseq(
    features=
        file.path(qza.16s.path, "core_metrics16s/rarefied_table.qza"),
    tree="SI_pipeline/merged/16s/rooted-tree16s.qza",
    "SI_pipeline/merged/16s/taxonomy16s.qza",
    metadata = "SI_pipeline/merged/16s/maps/combined-map-2018-2021.txt"
)



physeq16sR0
## plot(physeq16sR0@phy_tree, show.tip.label = FALSE)

feature.2.tax.16s <-
    read.table("SI_pipeline/merged/16s/taxonomy.tsv", sep="\t",
               header=TRUE)

feature.2.tax.16s$Taxon <- paste("16s", feature.2.tax.16s$Taxon, sep=':')

## convert to a phylo class which is more useful downstream
tree.16sR0 <- phy_tree(physeq16sR0, errorIfNULL=TRUE)

## match the tip labs to the table with feature ID and Taxon
tree.16sR0$tip.label  <-  feature.2.tax.16s$Taxon[match(tree.16sR0$tip.label,
                                           feature.2.tax.16s$Feature.ID)]


## 10-24-2023 Rebecca is dropping all sequences that are only resolved to the first level D_0__Bacteria

# Identify tips with labels exactly matching '16s:D_0__Bacteria'
matching_tips <- grep('^16s:d__Bacteria$', tree.16sR0$tip.label)

# Drop the matching tips
tree.16sR0 <- drop.tip(tree.16sR0, matching_tips)

# # Drop the tips that are NA
tree.16sR0 <- drop.tip(tree.16sR0, tree.16sR0$tip.label[is.na(tree.16sR0$tip.label)])

#drop unassigned tips
unassigned_tips <- grep('^16s:Unassigned$', tree.16sR0$tip.label)

# Drop the tips that are NA
tree.16sR0 <- drop.tip(tree.16sR0, tree.16sR0$tip.label[is.na(tree.16sR0$tip.label)])

plot(tree.16sR0, show.tip.label = FALSE)
save(physeq16sR0, tree.16sR0, file= "../skyIslands/data/physeq16s.Rdata")
## ***********************************************************************
## 16s networks
## ***********************************************************************

finalASVtable <- read.csv(file.path(qza.16s.path, "final_asv_table/16s_final_asv_table.csv"), header=TRUE)

rownames(finalASVtable) <- finalASVtable[,1]

finalASVtable[,1] <- NULL

## drop the barcode cols
finalASVtable <- finalASVtable %>%
  dplyr::select(starts_with("d__"))

#fixing naming inconsistency
colnames(finalASVtable) <- paste("16s", colnames(finalASVtable), sep=':')

colnames(finalASVtable) <- gsub("\\.__", "", colnames(finalASVtable))

colnames(finalASVtable) <- gsub("\\.", "; ", colnames(finalASVtable))

indiv.comm.16sR0 <-
    bipartite::empty(catchDups(makeComm(taxonomy16sR0,
                                        feature.2.tax.16s)))
## saving this out so we can make the tree visualizations
save(indiv.comm.16sR0, file="../skyIslands/data/indiv.comm16sR0.Rdata")

finalASVtable <- finalASVtable/rowSums(finalASVtable)

finalASVtable <- as.matrix(finalASVtable)

bees.16s <- c(rownames(finalASVtable))

comms <- list(finalASVtable)

species.16s <- unique(unlist(sapply(comms, colnames)))

merged.comm.16s <- plyr::rbind.fill(lapply(comms, as.data.frame))

## check with number of columns against number of unique species
dim(merged.comm.16s)
length(species.16s)

rownames(merged.comm.16s) <- bees.16s


length(tree.16sR0$tip.label[tree.16sR0$tip.label %in% colnames(finalASVtable)])
labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(finalASVtable))]
print(length(labels_mismatch)) #85


# fixing subgroup 1 issue
tree.16sR0$tip.label[grep("Subgroup_1", tree.16sR0$tip.label)]
colnames(finalASVtable)[grep("Subgroup_1", colnames(finalASVtable))]

colnames(finalASVtable) <- gsub("f__Acidobacteriaceae_; Subgroup_1; ", "f__Acidobacteriaceae_(Subgroup_1)", colnames(finalASVtable))
colnames(finalASVtable)[grep("Subgroup_1", colnames(finalASVtable))]

# check labels again
labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(finalASVtable))]
print(length(labels_mismatch)) #80

# fixing apibacter_sp. issue
tree.16sR0$tip.label[grep("Apibacter_sp.", tree.16sR0$tip.label)]
colnames(finalASVtable)[grep("Apibacter_sp", colnames(finalASVtable))]

colnames(finalASVtable) <- gsub("_sp; ", "_sp.", colnames(finalASVtable))

# check labels again
labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(finalASVtable))]
print(length(labels_mismatch)) #67

# fixing actinobacter issue
tree.16sR0$tip.label[grep("uncultured_actinobacterium", tree.16sR0$tip.label)]
colnames(finalASVtable)[grep("uncultured_actinobacterium", colnames(finalASVtable))]


colnames(finalASVtable) <- gsub("c__Actinobacteria; o__0319; 7L14; f__0319; 7L14; g__0319; 7L14; ", "c__Actinobacteria; o__0319-7L14; f__0319-7L14; g__0319-7L14; ", colnames(finalASVtable))
# check labels again
labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(finalASVtable))]
print(length(labels_mismatch)) #65

# fixing marine group issue
tree.16sR0$tip.label[grep("marine_group", tree.16sR0$tip.label)]
colnames(finalASVtable)[grep("marine_group", colnames(finalASVtable))]
colnames(finalASVtable) <- gsub("g__CL500; 29_marine_group", "g__CL500-29_marine_group", colnames(finalASVtable))

# check labels again
labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(finalASVtable))]
print(length(labels_mismatch)) #64

# fixing clostridia issue
tree.16sR0$tip.label[grep("g__W5053", tree.16sR0$tip.label)]
colnames(finalASVtable)[grep("g__W5053", colnames(finalASVtable))]
colnames(finalASVtable) <- gsub("o__Peptostreptococcales; Tissierellales; f__Peptostreptococcales; Tissierellales",
                                "o__Peptostreptococcales-Tissierellales; f__Peptostreptococcales-Tissierellales",
                                colnames(finalASVtable))

# check labels again
labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(finalASVtable))]
print(length(labels_mismatch)) #62

# fixing hafnia issue
tree.16sR0$tip.label[grep("Obesumbacterium", tree.16sR0$tip.label)]
colnames(finalASVtable)[grep("Obesumbacterium", colnames(finalASVtable))]
colnames(finalASVtable) <- gsub("g__Hafnia; Obesumbacterium",
                                "g__Hafnia-Obesumbacterium",
                                colnames(finalASVtable))

# check labels again
labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(finalASVtable))]
print(length(labels_mismatch)) #57

# fixing methylobacter issue
tree.16sR0$tip.label[grep("cerastii", tree.16sR0$tip.label)]
colnames(finalASVtable)[grep("cerastii", colnames(finalASVtable))]
colnames(finalASVtable) <- gsub("g__Methylobacterium; Methylorubrum",
                                "g__Methylobacterium-Methylorubrum",
                                colnames(finalASVtable))

# check labels again
labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(finalASVtable))]
print(length(labels_mismatch)) #37

# fixing uncultured soil issue
tree.16sR0$tip.label[grep("uncultured_soil", tree.16sR0$tip.label)]
colnames(finalASVtable)[grep("uncultured_soil", colnames(finalASVtable))]
colnames(finalASVtable) <- gsub("g__1174; 901; 12",
                                "g__1174-901-12",
                                colnames(finalASVtable))

# check labels again
labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(finalASVtable))]
print(length(labels_mismatch)) #27

# fixing allorhizobium issue
tree.16sR0$tip.label[grep("g__Allorhizobium", tree.16sR0$tip.label)]
colnames(finalASVtable)[grep("g__Allorhizobium", colnames(finalASVtable))]
colnames(finalASVtable) <- gsub("g__Allorhizobium; Neorhizobium; Pararhizobium; Rhizobium",
                                "g__Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
                                colnames(finalASVtable))

# check labels again
labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(finalASVtable))]
print(length(labels_mismatch)) #12

# fixing rhizobiales issue
tree.16sR0$tip.label[grep("Rhizobiales_bacterium", tree.16sR0$tip.label)]
colnames(finalASVtable)[grep("Rhizobiales_bacterium", colnames(finalASVtable))]
colnames(finalASVtable) <- gsub("f__D05; 2; g__D05; 2",
                                "f__D05-2; g__D05-2",
                                colnames(finalASVtable))

# check labels again
labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(finalASVtable))]
print(length(labels_mismatch)) #11

# fixing burkholderia issue
tree.16sR0$tip.label[grep("uncultured_Burkholderia", tree.16sR0$tip.label)]
colnames(finalASVtable)[grep("uncultured_Burkholderia", colnames(finalASVtable))]
colnames(finalASVtable) <- gsub("f__TRA3; 20; g__TRA3; 20",
                                "f__TRA3-20; g__TRA3-20",
                                colnames(finalASVtable))
colnames(finalASVtable) <- gsub("f__SC; I; 84; g__SC; I; 84",
                                "f__SC-I-84; g__SC-I-84",
                                colnames(finalASVtable))


# check labels again
labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(finalASVtable))]
print(length(labels_mismatch)) #7

# fixing paraburkholderia issue
tree.16sR0$tip.label[grep("Paraburkholderia", tree.16sR0$tip.label)]
colnames(finalASVtable)[grep("Paraburkholderia", colnames(finalASVtable))]
colnames(finalASVtable) <- gsub("g__Burkholderia; Caballeronia; Paraburkholderia",
                                "g__Burkholderia-Caballeronia-Paraburkholderia",
                                colnames(finalASVtable))

# check labels again
labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(finalASVtable))]
print(length(labels_mismatch)) #1

# fixing numbering hyphen issue
tree.16sR0$tip.label[grep("o__0319", tree.16sR0$tip.label)]
colnames(finalASVtable)[grep("o__0319", colnames(finalASVtable))]
colnames(finalASVtable) <- gsub("o__0319; 6G20; f__0319; 6G20; g__0319; 6G20",
                                "o__0319-6G20; f__0319-6G20; g__0319-6G20",
                                colnames(finalASVtable))

# check labels again
labels_mismatch <- tree.16sR0$tip.label[!(tree.16sR0$tip.label %in% colnames(finalASVtable))]
print(length(labels_mismatch)) #0! yay

save(finalASVtable, file= "../skyIslands/data/presAbsTable.Rdata")



## ***********************************************************************
## working with a merged 16s tree
## ***********************************************************************

# upload our mega 16s phylogenetic tree

mega16sdata <- read_qza("SI_pipeline/merged/16s/rooted-tree16s.qza")

tree.16s <- mega16sdata$data

## function from: https://stackoverflow.com/questions/38570074/phylogenetics-in-r-collapsing-descendant-tips-of-an-internal-node

drop_dupes <- function(tree, thres=1e-5){
  tips <- which(tree$edge[,2] %in% 1:Ntip(tree))
  toDrop <- tree$edge.length[tips] < thres
  drop.tip(tree,tree$tip.label[toDrop])
}

## tree.16s <- tip_glom(tree.16s, h=0.05)

## what is a good cutoff?
tree.16s <- drop_dupes(tree.16s, thres=1e-3)

# make distance matrix

tree.16s$tip.label  <-  feature.2.tax.16s$Taxon[match(tree.16s$tip.label,
                                           feature.2.tax.16s$Feature.ID)]

tree.16s <- drop.tip(tree.16s, which(duplicated(tree.16s$tip.label)))



## ***********************************************************************
## RBCL
## ***********************************************************************
## R0
# weightedUFrbclqzaR0 <- read_qza(file.path(qza.rbcl.path,
#                                           'core_metricsRBCLR0/weighted_unifrac_distance_matrix.qza'))
# 
# unweightedUFrbclqzaR0 <- read_qza(file.path(qza.rbcl.path,
#                                             'core_metricsRBCLR0/unweighted_unifrac_distance_matrix.qza'))
# 
# ## access data inside artifacts
# 
# # RBCL
# # R0
# wphylo.dist.rbclR0 <- weightedUFrbclqzaR0$data
# phylo.dist.rbclR0 <- unweightedUFrbclqzaR0$data
# 
# 
# #Note, when taxonomy is imported, a single string is returned along
# #with a confidence score.  For many analysis we will want to break up
# #this string and for that purpose the parse_taxonomy() function is
# #provided:
# 
# # R0
# taxonomyRBCLR0 <-
#   read_qza("SI_pipeline/merged/RBCL/final/core_metricsRBCLR0/rarefied_table.qza")
# taxonomyRBCLR0 <- taxonomyRBCLR0$data
# 
# 
# ## in the future do the same for other runs
# ## taxonomyRBCLR1 <- read_qza("SI_pipeline/merged/RBCL/final/core_metricsRBCLR1/rarefied_table.qza")
# ## taxonomyRBCLR1 <- taxonomyRBCLR1$data
# 
# 
# #Phyloseq tutorials here: https://joey711.github.io/phyloseq/
# 
# # ****CHECK PATHS BELOW ONCE FINISH RBCL PIPELINE****
# 
# ## NEED TO CHECK whether this actually has all of the samples and tips
# ## TODO: I think I need to rerun the pipeline steps to make sure the correct qzas are used
# 
# #RBCL 2023 phylogeny
# physeqRBCLR0 <- qza_to_phyloseq(
#   features="SI_pipeline/R2023/merged/RBCL/final/core_metricsRBCLR0/rarefied_table.qza",
#   tree="SI_pipeline/R2023/merged/RBCL/rooted-treeRBCL.qza",
#   "SI_pipeline/R2023/merged/RBCL/taxonomyRBCL.qza",
#   metadata = "SI_pipeline/R2023/merged/RBCL/maps/sky2020mapRBCL_combined_repsremoved.txt"
# )
# physeqRBCLR0
# plot(physeqRBCLR0@phy_tree, show.tip.label = FALSE)
# 
# #RBCL R1 phylogeny
# physeqRBCLR1 <- qza_to_phyloseq(
#   features="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR1/rarefied_table.qza",
#   tree="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/rooted-treeRBCL.qza",
#   "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCL.qza",
#   metadata = "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/maps/sky2020mapRBCL_combined_repsremoved.txt"
# )
# physeqRBCLR1
# plot(physeqRBCLR1@phy_tree, show.tip.label = FALSE)
# 
# #RBCL R2 phylogeny
# physeqRBCLR2 <- qza_to_phyloseq(
#   features="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR2/rarefied_table.qza",
#   tree="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/rooted-treeRBCL.qza",
#   "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCL.qza",
#   metadata = "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/maps/sky2020mapRBCL_combined_repsremoved.txt"
# )
# physeqRBCLR2
# plot(physeqRBCLR2@phy_tree, show.tip.label = FALSE)
# 
# #RBCL R3 phylogeny
# physeqRBCLR3 <- qza_to_phyloseq(
#   features="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR3/rarefied_table.qza",
#   tree="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/rooted-treeRBCL.qza",
#   "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCL.qza",
#   metadata = "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/maps/sky2020mapRBCL_combined_repsremoved.txt"
# )
# physeqRBCLR3
# plot(physeqRBCLR3@phy_tree, show.tip.label = FALSE)
# 
# #RBCL R4 phylogeny
# physeqRBCLR4 <- qza_to_phyloseq(
#   features="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR4/rarefied_table.qza",
#   tree="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/rooted-treeRBCL.qza",
#   "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCL.qza",
#   metadata = "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/maps/sky2020mapRBCL_combined_repsremoved.txt"
# )
# physeqRBCLR4
# plot(physeqRBCLR4@phy_tree, show.tip.label = FALSE)
# 
# #RBCL R5 phylogeny
# physeqRBCLR5 <- qza_to_phyloseq(
#   features="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR5/rarefied_table.qza",
#   tree="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/rooted-treeRBCL.qza",
#   "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCL.qza",
#   metadata = "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/maps/sky2020mapRBCL_combined_repsremoved.txt"
# )
# physeqRBCLR5
# plot(physeqRBCLR5@phy_tree, show.tip.label = FALSE)
# 
# #RBCL R6 phylogeny
# physeqRBCLR6 <- qza_to_phyloseq(
#   features="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR6/rarefied_table.qza",
#   tree="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/rooted-treeRBCL.qza",
#   "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCL.qza",
#   metadata = "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/maps/sky2020mapRBCL_combined_repsremoved.txt"
# )
# physeqRBCLR6
# plot(physeqRBCLR6@phy_tree, show.tip.label = FALSE)
# 
# #RBCL R7 phylogeny
# physeqRBCLR7 <- qza_to_phyloseq(
#   features="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR7/rarefied_table.qza",
#   tree="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/rooted-treeRBCL.qza",
#   "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCL.qza",
#   metadata = "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/maps/sky2020mapRBCL_combined_repsremoved.txt"
# )
# physeqRBCLR7
# plot(physeqRBCLR7@phy_tree, show.tip.label = FALSE)
# ## in the future do the same for other runs
# ## #RBCL R1 phylogeny
# ## physeqRBCLR1 <- qza_to_phyloseq(
# ##   features="SI_pipeline/merged/RBCL/final/core_metricsRBCLR1/rarefied_table.qza",
# ##   tree="SI_pipeline/merged/RBCL/rooted-treeRBCL.qza",
# ##   "SI_pipeline/merged/RBCL/taxonomyRBCL.qza",
# ##   metadata = "SI_pipeline/merged/RBCL/maps/SI2019_R1mapRBCL.txt"
# ## )
# ## physeqRBCLR1
# 
# ## plot(physeqRBCLR1@phy_tree, show.tip.label = FALSE)
# ## ***********************************************************************
# 
# feature.2.tax.rbcl <-
#     read.table("SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCL.txt", sep="\t",
#                header=TRUE)
# 
# ## ## R knows to ignore the row in the txt with a # (see defaults read.table)
# 
# ## add 16s/rbcl to make the columns easy to find when they are added
# ## to the specimen data
# 
# feature.2.tax.rbcl$Taxon  <- gsub(" ", "_", feature.2.tax.rbcl$Taxon)
# feature.2.tax.rbcl$Taxon <- paste("RBCL", feature.2.tax.rbcl$Taxon,
#                                   sep=':')
# ## R0
# tree.rbclR0 <- phy_tree(physeqRBCLR0, errorIfNULL=TRUE)
# ## R1
# tree.rbclR1 <- phy_tree(physeqRBCLR1, errorIfNULL=TRUE)
# ## R2
# tree.rbclR2 <- phy_tree(physeqRBCLR2, errorIfNULL=TRUE)
# ## R3
# tree.rbclR3 <- phy_tree(physeqRBCLR3, errorIfNULL=TRUE)
# ## R4
# tree.rbclR4 <- phy_tree(physeqRBCLR4, errorIfNULL=TRUE)
# ## R5
# tree.rbclR5 <- phy_tree(physeqRBCLR5, errorIfNULL=TRUE)
# ## R6
# tree.rbclR6 <- phy_tree(physeqRBCLR6, errorIfNULL=TRUE)
# ## R7
# tree.rbclR7 <- phy_tree(physeqRBCLR7, errorIfNULL=TRUE)
# 
# ## in the future do the same for other runs
# ## tree.rbclR1 <- phy_tree(physeqRBCLR1, errorIfNULL=TRUE)
# 
# ## match the tip labs to the table with feature ID and Taxon
# ## R0
# tree.rbclR0$tip.label  <-  feature.2.tax.rbcl$Taxon[
#                            match(tree.rbclR0$tip.label,
#                            feature.2.tax.rbcl$Feature.ID)]
# ## R1
# tree.rbclR1$tip.label  <-  feature.2.tax.rbcl$Taxon[
#   match(tree.rbclR1$tip.label,
#         feature.2.tax.rbcl$Feature.ID)]
# ## R2
# tree.rbclR2$tip.label  <-  feature.2.tax.rbcl$Taxon[
#   match(tree.rbclR2$tip.label,
#         feature.2.tax.rbcl$Feature.ID)]
# ## R3
# tree.rbclR3$tip.label  <-  feature.2.tax.rbcl$Taxon[
#   match(tree.rbclR3$tip.label,
#         feature.2.tax.rbcl$Feature.ID)]
# ## R4
# tree.rbclR4$tip.label  <-  feature.2.tax.rbcl$Taxon[
#   match(tree.rbclR4$tip.label,
#         feature.2.tax.rbcl$Feature.ID)]
# ## R5
# tree.rbclR5$tip.label  <-  feature.2.tax.rbcl$Taxon[
#   match(tree.rbclR5$tip.label,
#         feature.2.tax.rbcl$Feature.ID)]
# ## R6
# tree.rbclR6$tip.label  <-  feature.2.tax.rbcl$Taxon[
#   match(tree.rbclR6$tip.label,
#         feature.2.tax.rbcl$Feature.ID)]
# ## R7
# tree.rbclR7$tip.label  <-  feature.2.tax.rbcl$Taxon[
#   match(tree.rbclR7$tip.label,
#         feature.2.tax.rbcl$Feature.ID)]
# 
# ## in the future do the same for other runs
# 
# ## tree.rbclR1$tip.label  <-  feature.2.tax.rbcl$Taxon[
# ##                            match(tree.rbclR1$tip.label,
# ##                            feature.2.tax.rbcl$Feature.ID)]
# 
# ## ***********************************************************************
# ## rbcl networks
# ## ***********************************************************************
# 
# #R0
# indiv.comm.rbclR0 <-
#     bipartite::empty(catchDups(makeComm(taxonomyRBCLR0,
#                                         feature.2.tax.rbcl,
#                                         feature.col="Feature.ID")))
#                                         # changed to "Feature.ID" to match 'makeComm' function
# 
# indiv.comm.rbclR0 <- indiv.comm.rbclR0/rowSums(indiv.comm.rbclR0)
# ## R1
# indiv.comm.rbclR1 <-
#   bipartite::empty(catchDups(makeComm(taxonomyRBCLR1,
#                                       feature.2.tax.rbcl,
#                                       feature.col="Feature.ID")))
# 
# indiv.comm.rbclR1 <- indiv.comm.rbclR1/rowSums(indiv.comm.rbclR1)
# ## R2
# indiv.comm.rbclR2 <-
#   bipartite::empty(catchDups(makeComm(taxonomyRBCLR2,
#                                       feature.2.tax.rbcl,
#                                       feature.col="Feature.ID")))
# 
# indiv.comm.rbclR2 <- indiv.comm.rbclR2/rowSums(indiv.comm.rbclR2)
# ## R3
# indiv.comm.rbclR3 <-
#   bipartite::empty(catchDups(makeComm(taxonomyRBCLR3,
#                                       feature.2.tax.rbcl,
#                                       feature.col="Feature.ID")))
# 
# indiv.comm.rbclR3 <- indiv.comm.rbclR3/rowSums(indiv.comm.rbclR3)
# ## R4
# indiv.comm.rbclR4 <-
#   bipartite::empty(catchDups(makeComm(taxonomyRBCLR4,
#                                       feature.2.tax.rbcl,
#                                       feature.col="Feature.ID")))
# 
# indiv.comm.rbclR4 <- indiv.comm.rbclR4/rowSums(indiv.comm.rbclR4)
# ## R5
# indiv.comm.rbclR5 <-
#   bipartite::empty(catchDups(makeComm(taxonomyRBCLR5,
#                                       feature.2.tax.rbcl,
#                                       feature.col="Feature.ID")))
# 
# indiv.comm.rbclR5 <- indiv.comm.rbclR5/rowSums(indiv.comm.rbclR5)
# ## R6
# indiv.comm.rbclR6 <-
#   bipartite::empty(catchDups(makeComm(taxonomyRBCLR6,
#                                       feature.2.tax.rbcl,
#                                       feature.col="Feature.ID")))
# 
# indiv.comm.rbclR6 <- indiv.comm.rbclR6/rowSums(indiv.comm.rbclR6)
# ## R7
# indiv.comm.rbclR7 <-
#   bipartite::empty(catchDups(makeComm(taxonomyRBCLR7,
#                                       feature.2.tax.rbcl,
#                                       feature.col="Feature.ID")))
# 
# indiv.comm.rbclR7 <- indiv.comm.rbclR7/rowSums(indiv.comm.rbclR7)
# 
# ## indiv.comm.rbclR1 <-
# ##     bipartite::empty(catchDups(makeComm(taxonomyRBCLR1,
# ##                                         feature.2.tax.rbcl)))
# ## indiv.comm.rbclR1 <- indiv.comm.rbclR1/rowSums(indiv.comm.rbclR1)
# 
# 
# 
# # bees.rbcl <- c(rownames(indiv.comm.rbclR0))
# 
# # comms <- list(indiv.comm.rbclR0)
# 
#  bees.rbcl <- c(rownames(indiv.comm.rbclR0),
#                rownames(indiv.comm.rbclR1),
#                rownames(indiv.comm.rbclR2),
#                rownames(indiv.comm.rbclR3),
#                rownames(indiv.comm.rbclR4),
#                rownames(indiv.comm.rbclR5),
#                rownames(indiv.comm.rbclR6),
#                rownames(indiv.comm.rbclR7))
# 
#  comms <- list(indiv.comm.rbclR0, indiv.comm.rbclR1,
#                indiv.comm.rbclR2, indiv.comm.rbclR3,
#                indiv.comm.rbclR4, indiv.comm.rbclR5,
#                indiv.comm.rbclR6, indiv.comm.rbclR7)
# 
# species.rbcl <- unique(unlist(sapply(comms, colnames)))
# 
# merged.comm.rbcl <- plyr::rbind.fill(lapply(comms, as.data.frame))
# 
# ## check with number of columns against number of unique species
# dim(merged.comm.rbcl)
# length(species.rbcl)
# 
# rownames(merged.comm.rbcl) <- bees.rbcl
# 
# megaRBCLdata <- read_qza("SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/rooted-treeRBCL.qza")
# tree.rbcl <- megaRBCLdata$data
# 
# ## tree.rbcl <- tip_glom(tree.rbcl, h=0.1)
# 
# ## tree.rbcl <- drop_dupes(tree.rbcl, thres=1e-5)
# 
# tree.rbcl$tip.label  <-  feature.2.tax.rbcl$Taxon[match(tree.rbcl$tip.label,
#                                            feature.2.tax.rbcl$Feature.ID)]


## ***********************************************************************
## Make mega dataset
## ***********************************************************************

## spec already includes parasite data
load('../skyIslands/data/spec_net_fdiv.Rdata')

#indiv.comm.rbcl <- as.data.frame(merged.comm.rbcl)
#pollen <- colnames(indiv.comm.rbcl)
#indiv.comm.rbcl$UniqueID <- rownames(indiv.comm.rbcl)

#16s
indiv.comm.16s <- as.data.frame(merged.comm.16s)
bact <- colnames(indiv.comm.16s)
indiv.comm.16s$UniqueID <- rownames(indiv.comm.16s)

## 4-9-2025 no duplicae UniqueIDs here

spec.net <-cbind(spec.net, indiv.comm.16s[, bact][match(spec.net$UniqueID,
                           indiv.comm.16s$UniqueID),])

## 4-9-2025 no duplicae UniqueIDs here

#spec.net <-cbind(spec.net, indiv.comm.rbcl[, pollen][match(spec.net$UniqueID,
#                           indiv.comm.rbcl$UniqueID),])

## adding in PD code to merge to spec.net
## FULL MICROBE DATASET

## pull out 16s columns
microbes <- colnames(spec.net)[grepl("16s:", colnames(spec.net))] 
microbes <- microbes[microbes %in% tree.16s$tip.label]

## check which only have NAs in these columns (not screened) and drop them
screened.microbes <- apply(spec.net, 1, function(x) all(is.na(x[microbes])))
spec.microbes <- spec.net[!screened.microbes, ]

## 4-9-2025 no duplicae UniqueIDs here

#3 rows have 0 for all microbes, need to drop
spec.microbes <- spec.microbes[rowSums(spec.microbes[,microbes])!=0,]
## 4-9-2025 no duplicae UniqueIDs here

## QUESTION: should include root = TRUE? if false gives warning 3x
## warning: Rooted tree and include.root=TRUE argument required to calculate PD of single-species communities. Single species community assigned PD value of NA.

## phylogenetic distance function, modified from picante
PD <- apply(spec.microbes[,microbes], 1, function(x){
  this.bee <- x[x > 0]
  this.tree <- picante::prune.sample(t(this.bee), tree.16s)
  picante::pd(t(this.bee), this.tree, include.root = TRUE)
  #browser()
})

## phylogenetic species evenness for each bee
PSE <- apply(spec.microbes[,microbes], 1, function(x){
  this.bee <- x[x > 0]
  this.tree <- picante::prune.sample(t(this.bee), tree.16s)
  picante::pse(t(this.bee), this.tree)
  #browser()
})

## psd() phylogenetic species diversity metrics/phylogenetic species richness (PSR)
## ses.pd standardized effect size of phylo diversity in communities 

PD <- do.call(rbind, PD)
spec.microbes <- cbind(spec.microbes, PD)

## 4-9-2025 no duplicae UniqueIDs here

PSE <- do.call(rbind, PSE) %>%
  select(PSEs)

spec.microbes <- cbind(spec.microbes, PSE)
## 4-9-2025 no duplicae UniqueIDs here

## Merge back onto specimen data
spec.net <- merge(spec.net, spec.microbes, all.x=TRUE)
## 4-9-2025 no duplicae UniqueIDs here
spec.net$ScreenedMicrobes <- ifelse(spec.net$Site %in% unique(spec.microbes$Site), 1, 0)
## 4-9-2025 no duplicae UniqueIDs here
## change microbe NAs to 0
spec.net <- spec.net %>%
  mutate(PD = replace_na(PD, 0),
         PSEs = replace_na(PSEs, 0))
## 4-9-2025 no duplicae UniqueIDs here
spec.net[,microbes][is.na(spec.net[,microbes])] <- 0
## 4-9-2025 no duplicae UniqueIDs here

## ONLY OBLIGATE MICROBES DATASET


## splitting out obligate bee microbes based on Zheng and Moran paper
bee.obligates <- "Lactobacillaceae|Bifidobacteriaceae|Neisseriaceae|Orbaceae|Bartonella|Acetobacteraceae"

## this is a list of the microbe strains that contain the known obligate bee microbe genera
bee.obligate.microbes <- microbes[grepl(bee.obligates, microbes, fixed=FALSE)]
bee.obligate.microbes <- bee.obligate.microbes[bee.obligate.microbes %in% tree.16s$tip.label]
## now need to subset spec.microbes to be just the microbe columns with bee obligates and calculate PD

## phylogenetic distance function, modified from picante
PD.obligate <- apply(spec.microbes[,bee.obligate.microbes], 1, function(x){
  tryCatch({
    this.bee <- x[x > 0]
    this.tree <- picante::prune.sample(t(this.bee), tree.16s)
    pd_value <- picante::pd(t(this.bee), this.tree, include.root = TRUE)
    if (is.null(pd_value) || length(pd_value) == 0) pd_value <- 0  # Assign zero if PD is NULL or empty list
    else pd_value[[1]]  # Extract the first element if PD is a list
    data.frame(PD = pd_value[[1]], SR = pd_value[[2]])
  }, error = function(e) {
    if (grepl("Tree has no branch lengths", e$message)) {
      # If error is due to tree having no branch lengths, return zero for PD and SR
      return(data.frame(PD = 0, SR = 0))
    } else {
      # If it's a different error, re-raise it
      stop(e)
    }
  })
})

# Convert the result into a dataframe
result_df <- do.call(rbind, PD.obligate) 

# Rename column names to indicate these are the obligate only pd and sr
names(result_df)[names(result_df) == "PD"] <- "PD.obligate"
names(result_df)[names(result_df) == "SR"] <- "SR.obligate"


spec.microbes <- cbind(spec.microbes, result_df) 


## Merge back onto specimen data
spec.net <- merge(spec.net, spec.microbes, all.x=TRUE)
## 4-9-2025 no duplicae UniqueIDs here
## change microbe NAs to 0
spec.net <- spec.net %>%
  mutate(PD.obligate = replace_na(PD.obligate, 0))

spec.net[,microbes][is.na(spec.net[,microbes])] <- 0
## 4-9-2025 no duplicae UniqueIDs here
## ONLY TRANSIENT MICROBES DATASET

## splitting out obligate bee microbes based on Zheng and Moran paper
bee.obligates <- "Lactobacillaceae|Bifidobacteriaceae|Neisseriaceae|Orbaceae|Bartonella|Acetobacteraceae"

## this is a list of the microbe strains that contain the known transient bee microbe genera
bee.transient.microbes <- microbes[!grepl(bee.obligates, microbes, fixed=FALSE)]
bee.transient.microbes <- bee.transient.microbes[bee.transient.microbes %in% tree.16s$tip.label]
## now need to subset spec.microbes to be just the microbe columns with bee transients and calculate PD

## phylogenetic distance function, modified from picante
PD.transient <- apply(spec.microbes[,bee.transient.microbes], 1, function(x){
  tryCatch({
    this.bee <- x[x > 0]
    this.tree <- picante::prune.sample(t(this.bee), tree.16s)
    pd_value <- picante::pd(t(this.bee), this.tree, include.root = TRUE)
    if (is.null(pd_value) || length(pd_value) == 0) pd_value <- 0  # Assign zero if PD is NULL or empty list
    else pd_value[[1]]  # Extract the first element if PD is a list
    data.frame(PD = pd_value[[1]], SR = pd_value[[2]])
  }, error = function(e) {
    if (grepl("Tree has no branch lengths", e$message)) {
      # If error is due to tree having no branch lengths, return zero for PD and SR
      return(data.frame(PD = 0, SR = 0))
    } else {
      # If it's a different error, re-raise it
      stop(e)
    }
  })
})

# Convert the result into a dataframe
trans_df <- do.call(rbind, PD.transient) 

# Rename column names to indicate these are the transient only pd and sr
names(trans_df)[names(trans_df) == "PD"] <- "PD.transient"
names(trans_df)[names(trans_df) == "SR"] <- "SR.transient"


spec.microbes <- cbind(spec.microbes, trans_df)
## 4-9-2025 no duplicae UniqueIDs here
## Merge back onto specimen data
spec.net <- merge(spec.net, spec.microbes, all.x=TRUE)
## 4-9-2025 no duplicae UniqueIDs here
## change microbe NAs to 0
spec.net <- spec.net %>%
  mutate(PD.transient = replace_na(PD.transient, 0))

spec.net[,microbes][is.na(spec.net[,microbes])] <- 0


## drop duplicate sampleID. LP: This should not happen?
## spec.net <- spec.net[!duplicated(spec.net$UniqueID),]
any(duplicated(spec.net$UniqueID))

save(spec.net, file= "../skyIslands/data/spec_RBCL_16s.Rdata")

write.csv(spec.net, file= "../skyIslands/data/spec_RBCL_16s.csv",
          row.names=FALSE)

save(tree.16s, #tree.rbcl,
     file="../skyIslands/data/trees.Rdata")


