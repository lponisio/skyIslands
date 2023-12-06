library(ggplot2)
library(brms)
library(dplyr)
library(tidyr)
library(tidybayes)
library(ggthemes)
library(bayestestR)


load('../../data/spec_RBCL_16s.Rdata')
load("../../data/trees.Rdata")
site.sum <- read.csv("../../data/sitestats.csv")


spec.net <- spec.net[!is.na(spec.net$GenusSpecies),]

spec.net <- spec.net[spec.net$Family != "Syrphidae",]

parasites <- c(#"AspergillusSpp", ## problematic parasite!
  "AscosphaeraSpp",
  "ApicystisSpp",
  "CrithidiaExpoeki",
  "CrithidiaBombi",
  "CrithidiaSpp",
  "NosemaBombi",
  "NosemaCeranae")

spec.net <- merge(spec.net, site.sum, all.x=TRUE)

traits <-
    read.csv("../../../skyIslands_saved/data/raw/bee_traits.csv")
traits$GenusSpecies <- fix.white.space(traits$GenusSpecies)
traits <- traits[, c("GenusSpecies", "Sociality", "Lecty", "MeanITD"),]

net.traits <- read.csv("../../data/networks_traits.csv")
net.traits <- net.traits[, c("GenusSpecies", "r.degree"),]

traits <- merge(traits, net.traits, by="GenusSpecies", all.x=TRUE)

spec.net <- merge(spec.net, traits, all.x=TRUE, by="GenusSpecies")



##load tree from :
##Henriquez Piskulich, Patricia Andrea; Hugall, Andrew F.; Stuart-Fox, Devi (2023). A supermatrix phylogeny of the world’s bees (Hymenoptera: Anthophila) [Dataset]. Dryad. https://doi.org/10.5061/dryad.80gb5mkw1

phylo <- ape::read.tree("../../data/BEE_mat7_fulltree.nwk")


## i think we need to drop the tips that aren't in our dataset
## changing sp labels to match phylogeny
spec.net$GenusSpecies <- gsub("Megachile gemula fulvogemula", "Megachile gemula", spec.net$GenusSpecies)
spec.net$GenusSpecies <- gsub("Megachile melanophaea rohweri", "Megachile melanophaea", spec.net$GenusSpecies)


##clean up unwanted portion of labels
pattern <- "(_n\\d+m\\d+_[A-Za-z0-9]+)?$"
phylo$tip.label <- gsub(pattern, "", phylo$tip.label)

## replace underscore with space
phylo$tip.label <- gsub("_", " ", phylo$tip.label)

## i think we need to drop the tips that aren't in our dataset
## changing sp labels to match phylogeny
species_to_keep <- na.omit(unique(spec.net$GenusSpecies))

## megachile comate and megachile subexilis are not in phylogeny so will drop these
species_to_keep <- species_to_keep[!species_to_keep %in% c("Megachile comata", "Megachile subexilis")]

phylo_tips <- phylo$tip.label
#only keep tips that match our species
phylo <- ape::keep.tip(phylo, species_to_keep[species_to_keep %in% phylo_tips])

phylo_matrix <- ape::vcv.phylo(phylo)

## dropping species not in the phylogeny from the dataset
drop.species <- unique(spec.net$GenusSpecies[!(spec.net$GenusSpecies %in% rownames(phylo_matrix))])
spec.net <- spec.net[!spec.net$GenusSpecies %in% drop.species,]

# # changing missing sociality to Unknown
# spec.net$Sociality <- spec.net$Sociality %>%
#   replace_na("Unknown")
# 
# # changing missing MeanITD to 0
# spec.net$MeanITD <- spec.net$MeanITD %>%
#   replace_na(0)
# 
# # changing missing rare.degree to 0
# spec.net$rare.degree <- spec.net$rare.degree %>%
#   replace_na(0)



dir.create(path="saved", showWarnings = FALSE)
dir.create(path="saved/tables", showWarnings = FALSE)

spec.net <- spec.net[order(spec.net$Site),]

## create a dumby varaible "Weight" to deal with the data sets being at
## different levels to get around the issue of having to pass in one
## data set into brms
spec.net$YearSR <- paste(spec.net$Year, spec.net$SampleRound, sep=";")
spec.net$YearSRGenusSpecies <- paste(spec.net$YearSR, spec.net$GenusSpecies, sep=";")

## will need to modify when we have multiple years
spec.net <- makeDataMultiLevel(spec.net, "Site", "YearSR")

spec.net[, variables.to.log] <- log(spec.net[,variables.to.log])

##  center all of the x variables, need to use unique values to avoid
##  repetition by the number of specimens
spec.net <- standardizeVars(spec.net, vars_yearsr, "YearSR")

#spec.net <- standardizeVars(spec.net, vars_sp, "YearSRGenusSpecies")
## create a dumby varaible "WeightPar" for the parasite data. The
## original intention was to keep stan from dropping data for
## site-level models, but weight is 0 for parasite models.
spec.net <- prepParasiteWeights()


#genus_pd_fit <- function(spec.net, this_genus, num_iter){

microbes <- colnames(spec.net)[grepl("16s:", colnames(spec.net))] 

screened.microbes <- apply(spec.net, 1, function(x) all(is.na(x[microbes])))

spec.microbes <- spec.net[!screened.microbes, ]

##adding abundance weights column
abund_csv <- data.frame(read.csv("../../data/sp_year_site_round.csv"))

#join abundance csv
spec.net <- join(spec.net, abund_csv)

#multiply weightspar by abundance to get abundance weights
spec.net$WeightsAbund <- spec.net$WeightsPar * spec.net$AbundanceSYR
spec.net$LogWeightsAbund <- log(spec.net$WeightsAbund + 1)


## QUESTION: should include root = TRUE? if false gives warning 3x
## warning: Rooted tree and include.root=TRUE argument required to calculate PD of single-species communities. Single species community assigned PD value of NA.
PD <- apply(spec.microbes[,microbes], 1, function(x){
  this.bee <- x[x > 0]
  this.tree <- prune.sample(t(this.bee), tree.16s)
  #browser()
  picante::pd(t(this.bee), this.tree, include.root = TRUE)
})

PD <- do.call(rbind, PD)

spec.microbes <- cbind(spec.microbes, PD)

spec.net <- merge(spec.net, spec.microbes, all.x=TRUE)

## check which individuals don't have microbe data
drop.PD.NA <- unique(spec.net$UniqueID[spec.net$WeightsPar == 1 &
                                         is.na(spec.net$PD)])

## drop individuals that had parasite screen but not microbe
## filling in zeros for NAs in PD
spec.net <- spec.net[!(spec.net$UniqueID %in% drop.PD.NA),] %>%
  mutate(PD = ifelse(!is.na(PD), PD, 0))


#genus.microbes <- spec.microbes[spec.microbes$Genus == this_genus, ]

spec.all <- spec.net


## bombus only data
spec.bombus <- spec.all
spec.bombus$LogWeightsAbund[spec.bombus$Genus != "Bombus"] <- 0

# 
# ## megachile comate and megachile subexilis are not in phylogeny so will drop these
# species_to_keep <- species_to_keep[!species_to_keep %in% c("Megachile comata", "Megachile subexilis")]
# 
# phylo_tips <- phylo$tip.label
# #only keep tips that match our species
# phylo <- ape::keep.tip(phylo, species_to_keep[species_to_keep %in% phylo_tips])
# 
# phylo_matrix <- ape::vcv.phylo(phylo)


## apis only data
spec.apis <- spec.all
spec.apis$LogWeightsAbund[spec.apis$Genus != "Apis"] <- 0

## melissodes only data
spec.melissodes <- spec.all
spec.melissodes$WeightsPar[spec.melissodes$Genus != "Melissodes"] <- 0

spec.apidae <- spec.all
spec.apidae$WeightsPar[spec.melissodes$Family != "Apidae"] <- 0


