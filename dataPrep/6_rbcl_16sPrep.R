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
setwd('skyIslands')
setwd('../skyIslands_saved')


library(tidyr)
library(dplyr)
library(bipartite)
library(phyloseq)
library(TreeTools)
library(devtools)

devtools::install_github("jbisanz/qiime2R")
library(qiime2R)

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

save(physeq16sR0, file= "../skyIslands/data/physeq16s.Rdata")

#physeq16sR0
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
# Load the required package
library(ape)


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




## ***********************************************************************
## 16s networks
## ***********************************************************************

## something going wrong here, the make indiv.comm.16sR0 object is coming back blank
## filtering step maybe not working?


#2018 samples 
comm.matrix <- makeComm(taxonomy16sR0,
                             feature.2.tax.16s)


deduplicated.matrix <- catchDups(comm.matrix)

indiv.comm.16sR0 <-
    bipartite::empty(catchDups(makeComm(taxonomy16sR0,
                                        feature.2.tax.16s)))
indiv.comm.16sR0 <- indiv.comm.16sR0/rowSums(indiv.comm.16sR0)

save(indiv.comm.16sR0, file= "../skyIslands/data/presAbsTable.Rdata")

bees.16s <- c(rownames(indiv.comm.16sR0))

comms <- list(indiv.comm.16sR0)

species.16s <- unique(unlist(sapply(comms, colnames)))

merged.comm.16s <- plyr::rbind.fill(lapply(comms, as.data.frame))

## check with number of columns against number of unique species
dim(merged.comm.16s)
length(species.16s)

rownames(merged.comm.16s) <- bees.16s



## ***********************************************************************
## working with a merged 16s tree
## ***********************************************************************
#2018 samples!!

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
# R0
weightedUFrbclqzaR0 <- read_qza(file.path(qza.rbcl.path,
                                          'core_metricsRBCLR0/weighted_unifrac_distance_matrix.qza'))

unweightedUFrbclqzaR0 <- read_qza(file.path(qza.rbcl.path,
                                            'core_metricsRBCLR0/unweighted_unifrac_distance_matrix.qza'))
# R1
weightedUFrbclqzaR1 <- read_qza(file.path(qza.rbcl.path,
                                          'core_metricsRBCLR1/weighted_unifrac_distance_matrix.qza'))

unweightedUFrbclqzaR1 <- read_qza(file.path(qza.rbcl.path,
                                            'core_metricsRBCLR1/unweighted_unifrac_distance_matrix.qza'))
# R2
weightedUFrbclqzaR2 <- read_qza(file.path(qza.rbcl.path,
                                          'core_metricsRBCLR2/weighted_unifrac_distance_matrix.qza'))

unweightedUFrbclqzaR2 <- read_qza(file.path(qza.rbcl.path,
                                            'core_metricsRBCLR2/unweighted_unifrac_distance_matrix.qza'))

# R3
weightedUFrbclqzaR3 <- read_qza(file.path(qza.rbcl.path,
                                          'core_metricsRBCLR3/weighted_unifrac_distance_matrix.qza'))

unweightedUFrbclqzaR3 <- read_qza(file.path(qza.rbcl.path,
                                            'core_metricsRBCLR3/unweighted_unifrac_distance_matrix.qza'))
# R4
weightedUFrbclqzaR4 <- read_qza(file.path(qza.rbcl.path,
                                          'core_metricsRBCLR4/weighted_unifrac_distance_matrix.qza'))

unweightedUFrbclqzaR4 <- read_qza(file.path(qza.rbcl.path,
                                            'core_metricsRBCLR4/unweighted_unifrac_distance_matrix.qza'))
# R5
weightedUFrbclqzaR5 <- read_qza(file.path(qza.rbcl.path,
                                          'core_metricsRBCLR5/weighted_unifrac_distance_matrix.qza'))

unweightedUFrbclqzaR5 <- read_qza(file.path(qza.rbcl.path,
                                            'core_metricsRBCLR5/unweighted_unifrac_distance_matrix.qza'))
# R6
weightedUFrbclqzaR6 <- read_qza(file.path(qza.rbcl.path,
                                          'core_metricsRBCLR6/weighted_unifrac_distance_matrix.qza'))

unweightedUFrbclqzaR6 <- read_qza(file.path(qza.rbcl.path,
                                            'core_metricsRBCLR6/unweighted_unifrac_distance_matrix.qza'))

# R7
weightedUFrbclqzaR7 <- read_qza(file.path(qza.rbcl.path,
                                          'core_metricsRBCLR7/weighted_unifrac_distance_matrix.qza'))

unweightedUFrbclqzaR7 <- read_qza(file.path(qza.rbcl.path,
                                            'core_metricsRBCLR7/unweighted_unifrac_distance_matrix.qza'))

## access data inside artifacts

# RBCL
# R0
wphylo.dist.rbclR0 <- weightedUFrbclqzaR0$data
phylo.dist.rbclR0 <- unweightedUFrbclqzaR0$data
# R1
wphylo.dist.rbclR1 <- weightedUFrbclqzaR1$data
phylo.dist.rbclR1 <- unweightedUFrbclqzaR1$data
# R2
wphylo.dist.rbclR2 <- weightedUFrbclqzaR2$data
phylo.dist.rbclR2 <- unweightedUFrbclqzaR2$data
# R3
wphylo.dist.rbclR3 <- weightedUFrbclqzaR3$data
phylo.dist.rbclR3 <- unweightedUFrbclqzaR3$data
# R4
wphylo.dist.rbclR4 <- weightedUFrbclqzaR4$data
phylo.dist.rbclR4 <- unweightedUFrbclqzaR4$data
# R5
wphylo.dist.rbclR5 <- weightedUFrbclqzaR5$data
phylo.dist.rbclR5 <- unweightedUFrbclqzaR5$data
# R6
wphylo.dist.rbclR6 <- weightedUFrbclqzaR6$data
phylo.dist.rbclR6 <- unweightedUFrbclqzaR6$data
# R7
wphylo.dist.rbclR7 <- weightedUFrbclqzaR7$data
phylo.dist.rbclR7 <- unweightedUFrbclqzaR7$data

## in the future do the same for other runs
## wphylo.dist.rbclR1 <- weightedUFrbclqzaR1$data
## phylo.dist.rbclR1 <- unweightedUFrbclqzaR1$data

#Note, when taxonomy is imported, a single string is returned along
#with a confidence score.  For many analysis we will want to break up
#this string and for that purpose the parse_taxonomy() function is
#provided:

# R0
taxonomyRBCLR0 <-
  read_qza("SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR0/rarefied_table.qza")
taxonomyRBCLR0 <- taxonomyRBCLR0$data

# R1
taxonomyRBCLR1 <-
  read_qza("SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR1/rarefied_table.qza")
taxonomyRBCLR1 <- taxonomyRBCLR1$data

# R2
taxonomyRBCLR2 <-
  read_qza("SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR2/rarefied_table.qza")
taxonomyRBCLR2 <- taxonomyRBCLR2$data

# R3
taxonomyRBCLR3 <-
  read_qza("SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR3/rarefied_table.qza")
taxonomyRBCLR3 <- taxonomyRBCLR3$data

# R4
taxonomyRBCLR4 <-
  read_qza("SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR4/rarefied_table.qza")
taxonomyRBCLR4 <- taxonomyRBCLR4$data

# R5
taxonomyRBCLR5 <-
  read_qza("SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR5/rarefied_table.qza")
taxonomyRBCLR5 <- taxonomyRBCLR5$data

# R6
taxonomyRBCLR6 <-
  read_qza("SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR6/rarefied_table.qza")
taxonomyRBCLR6 <- taxonomyRBCLR6$data

# R7
taxonomyRBCLR7 <-
  read_qza("SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR7/rarefied_table.qza")
taxonomyRBCLR7 <- taxonomyRBCLR7$data

## in the future do the same for other runs
## taxonomyRBCLR1 <- read_qza("SI_pipeline/merged/RBCL/final/core_metricsRBCLR1/rarefied_table.qza")
## taxonomyRBCLR1 <- taxonomyRBCLR1$data


#Phyloseq tutorials here: https://joey711.github.io/phyloseq/

# ****CHECK PATHS BELOW ONCE FINISH RBCL PIPELINE****

#RBCL R0 phylogeny
physeqRBCLR0 <- qza_to_phyloseq(
  features="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR0/rarefied_table.qza",
  tree="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/rooted-treeRBCL.qza",
  "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCL.qza",
  metadata = "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/maps/sky2020mapRBCL_combined_repsremoved.txt"
)
physeqRBCLR0
plot(physeqRBCLR0@phy_tree, show.tip.label = FALSE)

#RBCL R1 phylogeny
physeqRBCLR1 <- qza_to_phyloseq(
  features="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR1/rarefied_table.qza",
  tree="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/rooted-treeRBCL.qza",
  "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCL.qza",
  metadata = "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/maps/sky2020mapRBCL_combined_repsremoved.txt"
)
physeqRBCLR1
plot(physeqRBCLR1@phy_tree, show.tip.label = FALSE)

#RBCL R2 phylogeny
physeqRBCLR2 <- qza_to_phyloseq(
  features="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR2/rarefied_table.qza",
  tree="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/rooted-treeRBCL.qza",
  "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCL.qza",
  metadata = "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/maps/sky2020mapRBCL_combined_repsremoved.txt"
)
physeqRBCLR2
plot(physeqRBCLR2@phy_tree, show.tip.label = FALSE)

#RBCL R3 phylogeny
physeqRBCLR3 <- qza_to_phyloseq(
  features="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR3/rarefied_table.qza",
  tree="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/rooted-treeRBCL.qza",
  "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCL.qza",
  metadata = "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/maps/sky2020mapRBCL_combined_repsremoved.txt"
)
physeqRBCLR3
plot(physeqRBCLR3@phy_tree, show.tip.label = FALSE)

#RBCL R4 phylogeny
physeqRBCLR4 <- qza_to_phyloseq(
  features="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR4/rarefied_table.qza",
  tree="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/rooted-treeRBCL.qza",
  "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCL.qza",
  metadata = "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/maps/sky2020mapRBCL_combined_repsremoved.txt"
)
physeqRBCLR4
plot(physeqRBCLR4@phy_tree, show.tip.label = FALSE)

#RBCL R5 phylogeny
physeqRBCLR5 <- qza_to_phyloseq(
  features="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR5/rarefied_table.qza",
  tree="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/rooted-treeRBCL.qza",
  "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCL.qza",
  metadata = "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/maps/sky2020mapRBCL_combined_repsremoved.txt"
)
physeqRBCLR5
plot(physeqRBCLR5@phy_tree, show.tip.label = FALSE)

#RBCL R6 phylogeny
physeqRBCLR6 <- qza_to_phyloseq(
  features="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR6/rarefied_table.qza",
  tree="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/rooted-treeRBCL.qza",
  "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCL.qza",
  metadata = "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/maps/sky2020mapRBCL_combined_repsremoved.txt"
)
physeqRBCLR6
plot(physeqRBCLR6@phy_tree, show.tip.label = FALSE)

#RBCL R7 phylogeny
physeqRBCLR7 <- qza_to_phyloseq(
  features="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/final/core_metricsRBCLR7/rarefied_table.qza",
  tree="SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/rooted-treeRBCL.qza",
  "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCL.qza",
  metadata = "SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/maps/sky2020mapRBCL_combined_repsremoved.txt"
)
physeqRBCLR7
plot(physeqRBCLR7@phy_tree, show.tip.label = FALSE)
## in the future do the same for other runs
## #RBCL R1 phylogeny
## physeqRBCLR1 <- qza_to_phyloseq(
##   features="SI_pipeline/merged/RBCL/final/core_metricsRBCLR1/rarefied_table.qza",
##   tree="SI_pipeline/merged/RBCL/rooted-treeRBCL.qza",
##   "SI_pipeline/merged/RBCL/taxonomyRBCL.qza",
##   metadata = "SI_pipeline/merged/RBCL/maps/SI2019_R1mapRBCL.txt"
## )
## physeqRBCLR1

## plot(physeqRBCLR1@phy_tree, show.tip.label = FALSE)
## ***********************************************************************

feature.2.tax.rbcl <-
    read.table("SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCL.txt", sep="\t",
               header=TRUE)

## ## R knows to ignore the row in the txt with a # (see defaults read.table)

## add 16s/rbcl to make the columns easy to find when they are added
## to the specimen data

feature.2.tax.rbcl$Taxon  <- gsub(" ", "_", feature.2.tax.rbcl$Taxon)
feature.2.tax.rbcl$Taxon <- paste("RBCL", feature.2.tax.rbcl$Taxon,
                                  sep=':')
#R0
tree.rbclR0 <- phy_tree(physeqRBCLR0, errorIfNULL=TRUE)
# R1
tree.rbclR1 <- phy_tree(physeqRBCLR1, errorIfNULL=TRUE)
# R2
tree.rbclR2 <- phy_tree(physeqRBCLR2, errorIfNULL=TRUE)
# R3
tree.rbclR3 <- phy_tree(physeqRBCLR3, errorIfNULL=TRUE)
# R4
tree.rbclR4 <- phy_tree(physeqRBCLR4, errorIfNULL=TRUE)
# R5
tree.rbclR5 <- phy_tree(physeqRBCLR5, errorIfNULL=TRUE)
# R6
tree.rbclR6 <- phy_tree(physeqRBCLR6, errorIfNULL=TRUE)
# R7
tree.rbclR7 <- phy_tree(physeqRBCLR7, errorIfNULL=TRUE)

## in the future do the same for other runs
## tree.rbclR1 <- phy_tree(physeqRBCLR1, errorIfNULL=TRUE)

## match the tip labs to the table with feature ID and Taxon
# R0
tree.rbclR0$tip.label  <-  feature.2.tax.rbcl$Taxon[
                           match(tree.rbclR0$tip.label,
                           feature.2.tax.rbcl$Feature.ID)]
# R1
tree.rbclR1$tip.label  <-  feature.2.tax.rbcl$Taxon[
  match(tree.rbclR1$tip.label,
        feature.2.tax.rbcl$Feature.ID)]
# R2
tree.rbclR2$tip.label  <-  feature.2.tax.rbcl$Taxon[
  match(tree.rbclR2$tip.label,
        feature.2.tax.rbcl$Feature.ID)]
# R3
tree.rbclR3$tip.label  <-  feature.2.tax.rbcl$Taxon[
  match(tree.rbclR3$tip.label,
        feature.2.tax.rbcl$Feature.ID)]
# R4
tree.rbclR4$tip.label  <-  feature.2.tax.rbcl$Taxon[
  match(tree.rbclR4$tip.label,
        feature.2.tax.rbcl$Feature.ID)]
# R5
tree.rbclR5$tip.label  <-  feature.2.tax.rbcl$Taxon[
  match(tree.rbclR5$tip.label,
        feature.2.tax.rbcl$Feature.ID)]
# R6
tree.rbclR6$tip.label  <-  feature.2.tax.rbcl$Taxon[
  match(tree.rbclR6$tip.label,
        feature.2.tax.rbcl$Feature.ID)]
# R7
tree.rbclR7$tip.label  <-  feature.2.tax.rbcl$Taxon[
  match(tree.rbclR7$tip.label,
        feature.2.tax.rbcl$Feature.ID)]

## in the future do the same for other runs

## tree.rbclR1$tip.label  <-  feature.2.tax.rbcl$Taxon[
##                            match(tree.rbclR1$tip.label,
##                            feature.2.tax.rbcl$Feature.ID)]

## ***********************************************************************
## rbcl networks
## ***********************************************************************

#R0
indiv.comm.rbclR0 <-
    bipartite::empty(catchDups(makeComm(taxonomyRBCLR0,
                                        feature.2.tax.rbcl,
                                        feature.col="Feature.ID"))) # changed to "Feature.ID" to match 'makeComm' function

indiv.comm.rbclR0 <- indiv.comm.rbclR0/rowSums(indiv.comm.rbclR0)
# R1
indiv.comm.rbclR1 <-
  bipartite::empty(catchDups(makeComm(taxonomyRBCLR1,
                                      feature.2.tax.rbcl,
                                      feature.col="Feature.ID")))

indiv.comm.rbclR1 <- indiv.comm.rbclR1/rowSums(indiv.comm.rbclR1)
# R2
indiv.comm.rbclR2 <-
  bipartite::empty(catchDups(makeComm(taxonomyRBCLR2,
                                      feature.2.tax.rbcl,
                                      feature.col="Feature.ID")))

indiv.comm.rbclR2 <- indiv.comm.rbclR2/rowSums(indiv.comm.rbclR2)
# R3
indiv.comm.rbclR3 <-
  bipartite::empty(catchDups(makeComm(taxonomyRBCLR3,
                                      feature.2.tax.rbcl,
                                      feature.col="Feature.ID")))

indiv.comm.rbclR3 <- indiv.comm.rbclR3/rowSums(indiv.comm.rbclR3)
# R4
indiv.comm.rbclR4 <-
  bipartite::empty(catchDups(makeComm(taxonomyRBCLR4,
                                      feature.2.tax.rbcl,
                                      feature.col="Feature.ID")))

indiv.comm.rbclR4 <- indiv.comm.rbclR4/rowSums(indiv.comm.rbclR4)
# R5
indiv.comm.rbclR5 <-
  bipartite::empty(catchDups(makeComm(taxonomyRBCLR5,
                                      feature.2.tax.rbcl,
                                      feature.col="Feature.ID")))

indiv.comm.rbclR5 <- indiv.comm.rbclR5/rowSums(indiv.comm.rbclR5)
# R6
indiv.comm.rbclR6 <-
  bipartite::empty(catchDups(makeComm(taxonomyRBCLR6,
                                      feature.2.tax.rbcl,
                                      feature.col="Feature.ID")))

indiv.comm.rbclR6 <- indiv.comm.rbclR6/rowSums(indiv.comm.rbclR6)
# R7
indiv.comm.rbclR7 <-
  bipartite::empty(catchDups(makeComm(taxonomyRBCLR7,
                                      feature.2.tax.rbcl,
                                      feature.col="Feature.ID")))

indiv.comm.rbclR7 <- indiv.comm.rbclR7/rowSums(indiv.comm.rbclR7)

## indiv.comm.rbclR1 <-
##     bipartite::empty(catchDups(makeComm(taxonomyRBCLR1,
##                                         feature.2.tax.rbcl)))
## indiv.comm.rbclR1 <- indiv.comm.rbclR1/rowSums(indiv.comm.rbclR1)



# bees.rbcl <- c(rownames(indiv.comm.rbclR0))

# comms <- list(indiv.comm.rbclR0)

 bees.rbcl <- c(rownames(indiv.comm.rbclR0),
               paste0("SF", rownames(indiv.comm.rbclR1)),
               paste0("SF", rownames(indiv.comm.rbclR2)),
               paste0("SF", rownames(indiv.comm.rbclR3)),
               paste0("SF", rownames(indiv.comm.rbclR4)),
               paste0("SF", rownames(indiv.comm.rbclR5)),
               paste0("SF", rownames(indiv.comm.rbclR6)),
               paste0("SF", rownames(indiv.comm.rbclR7)))

 comms <- list(indiv.comm.rbclR0, indiv.comm.rbclR1,
               indiv.comm.rbclR2, indiv.comm.rbclR3,
               indiv.comm.rbclR4, indiv.comm.rbclR5,
               indiv.comm.rbclR6, indiv.comm.rbclR7)

species.rbcl <- unique(unlist(sapply(comms, colnames)))

merged.comm.rbcl <- plyr::rbind.fill(lapply(comms, as.data.frame))

## check with number of columns against number of unique species
dim(merged.comm.rbcl)
length(species.rbcl)

rownames(merged.comm.rbcl) <- bees.rbcl


megaRBCLdata <- read_qza("SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/rooted-treeRBCL.qza")
tree.rbcl <- megaRBCLdata$data

## tree.rbcl <- tip_glom(tree.rbcl, h=0.1)

## tree.rbcl <- drop_dupes(tree.rbcl, thres=1e-5)

tree.rbcl$tip.label  <-  feature.2.tax.rbcl$Taxon[match(tree.rbcl$tip.label,
                                           feature.2.tax.rbcl$Feature.ID)]


## ***********************************************************************
## Make mega dataset
## ***********************************************************************


## spec already includes parasite data
load('../skyIslands/data/spec_net.Rdata')


indiv.comm.rbcl <- as.data.frame(merged.comm.rbcl)
pollen <- colnames(indiv.comm.rbcl)
indiv.comm.rbcl$UniqueID <- rownames(indiv.comm.rbcl)

#2018 16s
indiv.comm.16s <- as.data.frame(merged.comm.16s)
bact <- colnames(indiv.comm.16s)
indiv.comm.16s$UniqueID <- rownames(indiv.comm.16s)


spec.net <-cbind(spec.net, indiv.comm.16s[, bact][match(spec.net$UniqueID,
                           indiv.comm.16s$UniqueID),])

spec.net <-cbind(spec.net, indiv.comm.rbcl[, pollen][match(spec.net$UniqueID,
                           indiv.comm.rbcl$UniqueID),])


#drop duplicate sampleIDs
spec.net <- spec.net[!duplicated(spec.net$UniqueID),]
  
save(spec.net, file= "../skyIslands/data/spec_RBCL_16s.Rdata")

write.csv(spec.net, file= "../skyIslands/data/spec_RBCL_16s.csv",
          row.names=FALSE)

save(tree.16s, tree.rbcl, file="../skyIslands/data/trees.Rdata")


