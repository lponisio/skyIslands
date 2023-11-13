## from this
## tutorial:https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121

## tutorial also have some good pcoa plots and heatmaps, which arent
## on here
#dir.bombus <- '/Volumes/bombus/Dropbox (University of Oregon)'

dir.bombus <- '~/Dropbox (University of Oregon)'

setwd(dir.bombus)
setwd('skyIslands_saved')
rm(list=ls())

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
qza.16s.path  <- "SI_pipeline/R2018/2023_sequence_results_raw/merged/16s/final"
qza.rbcl.path  <- "SI_pipeline/merged/RBCL/final"

# 16s
weightedUF16sqzaR0 <- read_qza(file.path(qza.16s.path,
                      'core_metrics16s-combined/weighted_unifrac_distance_matrix.qza'))

unweightedUF16sqzaR0 <- read_qza(file.path(qza.16s.path,
                      'core_metrics16s-combined/unweighted_unifrac_distance_matrix.qza'))

# 16s
wphylo.dist.16sR0 <- weightedUF16sqzaR0$data
phylo.dist.16sR0 <-unweightedUF16sqzaR0$data

## Note, when taxonomy is imported, a single string is returned along
## with a confidence score.  For many analysis we will want to break
## up this string and for that purpose the parse_taxonomy() function
## is provided:

taxonomy16sR0 <- read_qza(
    file.path(qza.16s.path, "core_metrics16s-combined/rarefied_table.qza"))

taxonomy16sR0 <- taxonomy16sR0$data

## ## 16s R0 phylogeny

physeq16sR0 <- qza_to_phyloseq(
    features=
        file.path(qza.16s.path, "core_metrics16s-combined/rarefied_table.qza"),
    tree="SI_pipeline/R2018/2023_sequence_results_raw/merged/16s/rooted-tree16s-combined.qza",
    "SI_pipeline/R2018/2023_sequence_results_raw/merged/16s/combined-taxonomy.qza",
    metadata = "SI_pipeline/R2018/2023_sequence_results_raw/merged/16s/maps/combined-map-2018-2021.txt"
)

#physeq16sR0
## plot(physeq16sR0@phy_tree, show.tip.label = FALSE)

feature.2.tax.16s <-
    read.table("SI_pipeline/R2018/2023_sequence_results_raw/merged/16s/1efec88b-38a0-454b-a858-11c17590b86e/data/taxonomy.tsv", sep="\t",
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
matching_tips <- grep('^16s:D_0__Bacteria$', tree.16sR0$tip.label)

# Drop the matching tips
tree.16sR0 <- drop.tip(tree.16sR0, matching_tips)

plot(tree.16sR0, show.tip.label = FALSE)




## ***********************************************************************
## 16s networks
## ***********************************************************************


#2018 samples 

indiv.comm.16sR0 <-
    bipartite::empty(catchDups(makeComm(taxonomy16sR0,
                                        feature.2.tax.16s)))
indiv.comm.16sR0 <- indiv.comm.16sR0/rowSums(indiv.comm.16sR0)

bees.16s <- c(rownames(indiv.comm.16sR0))

comms <- list(indiv.comm.16sR0)

## bees.16s <- c(rownames(indiv.comm.16sR0),
##               paste0("SF",  rownames(indiv.comm.16sR1)),
##               paste0("SF", rownames(indiv.comm.16sR2)),
##                      paste0("SF", rownames(indiv.comm.16sR3)),
##               paste0("SF", rownames(indiv.comm.16sR4)))

## comms <- list(indiv.comm.16sR0, indiv.comm.16sR1,
##               indiv.comm.16sR2, indiv.comm.16sR3,
##               indiv.comm.16sR4)

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

mega16sdata <- read_qza("SI_pipeline/R2018/2023_sequence_results_raw/merged/16s/rooted-tree16s-combined.qza")

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

weightedUFrbclqzaR0 <- read_qza(file.path(qza.rbcl.path,
                  'core_metricsRBCLR0/weighted_unifrac_distance_matrix.qza'))

unweightedUFrbclqzaR0 <- read_qza(file.path(qza.rbcl.path,
                 'core_metricsRBCLR0/unweighted_unifrac_distance_matrix.qza'))

## access data inside artifacts

# RBCL
wphylo.dist.rbclR0 <- weightedUFrbclqzaR0$data
phylo.dist.rbclR0 <- unweightedUFrbclqzaR0$data

## in the future do the same for other runs
## wphylo.dist.rbclR1 <- weightedUFrbclqzaR1$data
## phylo.dist.rbclR1 <- unweightedUFrbclqzaR1$data

#Note, when taxonomy is imported, a single string is returned along
#with a confidence score.  For many analysis we will want to break up
#this string and for that purpose the parse_taxonomy() function is
#provided:

taxonomyRBCLR0 <-
    read_qza("SI_pipeline/merged/RBCL/final/core_metricsRBCLR0/rarefied_table.qza")

taxonomyRBCLR0 <- taxonomyRBCLR0$data

## in the future do the same for other runs
## taxonomyRBCLR1 <- read_qza("SI_pipeline/merged/RBCL/final/core_metricsRBCLR1/rarefied_table.qza")
## taxonomyRBCLR1 <- taxonomyRBCLR1$data


#Phyloseq tutorials here: https://joey711.github.io/phyloseq/

# ****CHECK PATHS BELOW ONCE FINISH RBCL PIPELINE****

#RBCL R0 phylogeny
physeqRBCLR0 <- qza_to_phyloseq(
  features="SI_pipeline/merged/RBCL/final/core_metricsRBCLR0/rarefied_table.qza",
  tree="SI_pipeline/merged/RBCL/rooted-treeRBCL.qza",
  "SI_pipeline/merged/RBCL/taxonomyRBCL.qza",
  metadata = "SI_pipeline/merged/RBCL/maps/SI2018mapRBCL.txt"
)
physeqRBCLR0
plot(physeqRBCLR0@phy_tree, show.tip.label = FALSE)

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
    read.table("SI_pipeline/merged/RBCL/taxonomyRBCL.txt", sep="\t",
               header=TRUE)

## ## R knows to ignore the row in the txt with a # (see defaults read.table)

## add 16s/rbcl to make the columns easy to find when they are added
## to the specimen data

feature.2.tax.rbcl$Taxon  <- gsub(" ", "_", feature.2.tax.rbcl$Taxon)
feature.2.tax.rbcl$Taxon <- paste("RBCL", feature.2.tax.rbcl$Taxon,
                                  sep=':')

tree.rbclR0 <- phy_tree(physeqRBCLR0, errorIfNULL=TRUE)

## in the future do the same for other runs
## tree.rbclR1 <- phy_tree(physeqRBCLR1, errorIfNULL=TRUE)

## match the tip labs to the table with feature ID and Taxon
tree.rbclR0$tip.label  <-  feature.2.tax.rbcl$Taxon[
                           match(tree.rbclR0$tip.label,
                           feature.2.tax.rbcl$Feature.ID)]

## in the future do the same for other runs

## tree.rbclR1$tip.label  <-  feature.2.tax.rbcl$Taxon[
##                            match(tree.rbclR1$tip.label,
##                            feature.2.tax.rbcl$Feature.ID)]

## ***********************************************************************
## rbcl networks
## ***********************************************************************

indiv.comm.rbclR0 <-
    bipartite::empty(catchDups(makeComm(taxonomyRBCLR0,
                                        feature.2.tax.rbcl,
                                        feature.col="FeatureID")))

indiv.comm.rbclR0 <- indiv.comm.rbclR0/rowSums(indiv.comm.rbclR0)

## indiv.comm.rbclR1 <-
##     bipartite::empty(catchDups(makeComm(taxonomyRBCLR1,
##                                         feature.2.tax.rbcl)))
## indiv.comm.rbclR1 <- indiv.comm.rbclR1/rowSums(indiv.comm.rbclR1)



bees.rbcl <- c(rownames(indiv.comm.rbclR0))

comms <- list(indiv.comm.rbclR0)

## bees.rbcl <- c(rownames(indiv.comm.rbclR0),
##               paste0("SF",  rownames(indiv.comm.rbclR1)),
##               paste0("SF", rownames(indiv.comm.rbclR2)),
##                      paste0("SF", rownames(indiv.comm.rbclR3)),
##               paste0("SF", rownames(indiv.comm.rbclR4)))

## comms <- list(indiv.comm.rbclR0, indiv.comm.rbclR1,
##               indiv.comm.rbclR2, indiv.comm.rbclR3,
##               indiv.comm.rbclR4)

species.rbcl <- unique(unlist(sapply(comms, colnames)))

merged.comm.rbcl <- plyr::rbind.fill(lapply(comms, as.data.frame))

## check with number of columns against number of unique species
dim(merged.comm.rbcl)
length(species.rbcl)

rownames(merged.comm.rbcl) <- bees.rbcl


megaRBCLdata <- read_qza("SI_pipeline/merged/RBCL/rooted-treeRBCL.qza")
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


