library(ggplot2)
library(brms)
library(dplyr)
library(tidyr)
library(tidybayes)
library(ggthemes)
library(bayestestR)
library(gtools)

## plotting
library(tidyr)
library(dplyr)
library(viridis)
library(tidybayes)
library(gridExtra)
library(grid)
#library(scales)
library(RColorBrewer)

library(rstantools)
library(performance)
library(bayestestR)
#library(see)

save.dir <- "saved/tables"
if(!dir.exists(save.dir)) {
  dir.create(save.dir, showWarnings = FALSE)
}

fig.dir <- "figures"
if(!dir.exists(fig.dir)) {
  dir.create(save.dir, showWarnings = FALSE)
}



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

dir.create(path="saved", showWarnings = FALSE)
dir.create(path="saved/tables", showWarnings = FALSE)

spec.net <- spec.net[order(spec.net$Site),]


## TODO trying Nicole's function for prep SEM data
## all of the variables that are explanatory variables and thus need
## to be centered


# maybe change ones with no ID to NoID instead of Na
# TODO check if i need this
spec.net$GenusSpecies <- if_else(spec.net$GenusSpecies=='', 'NoID', spec.net$GenusSpecies)

## raw, non standardized data for plotting
spec.orig <- prepDataSEM(spec.net, variables.to.log, 
                         standardize=FALSE)

## Make SEM weights and standardize data.
spec.net <- prepDataSEM(spec.net, variables.to.log, 
                        vars_yearsr = vars_yearsr, vars_sp = vars_sp, 
                        vars_yearsrsp = vars_yearsrsp,
                        vars_site=vars_site)



##load tree from :
##Henriquez Piskulich, Patricia Andrea; Hugall, Andrew F.; Stuart-Fox, Devi (2023). A supermatrix phylogeny of the world’s bees (Hymenoptera: Anthophila) [Dataset]. Dryad. https://doi.org/10.5061/dryad.80gb5mkw1

phylo <- ape::read.tree("../../data/BEE_mat7_fulltree.nwk")


## i think we need to drop the tips that aren't in our dataset
## changing sp labels to match phylogeny
# spec.net$GenusSpecies <- gsub("Megachile gemula fulvogemula", "Megachile gemula", spec.net$GenusSpecies)
# spec.net$GenusSpecies <- gsub("Megachile melanophaea rohweri", "Megachile melanophaea", spec.net$GenusSpecies)


##clean up unwanted portion of labels
pattern <- "(_n\\d+m\\d+_[A-Za-z0-9]+)?$"
phylo$tip.label <- gsub(pattern, "", phylo$tip.label)

## replace underscore with space
phylo$tip.label <- gsub("_", " ", phylo$tip.label)

## TODO check if this works
## Species that are not in the phylogeny are not used. brms is not allowing an incomplete
## phylogeny, to avoid the error we changed the species not present to one that is in the phylogeny. 
## We chose a species for which we did not do parasite screening and should not influence results.
not_in_phylo <- unique(spec.net$GenusSpecies[!spec.net$GenusSpecies %in% phylo$tip.label])
spec.net$GenusSpecies[spec.net$GenusSpecies %in% not_in_phylo]<- "Agapostemon angelicus"


# ## i think we need to drop the tips that aren't in our dataset
# ## changing sp labels to match phylogeny
# species_to_keep <- na.omit(unique(spec.net$GenusSpecies))
# 
# ## megachile comate and megachile subexilis are not in phylogeny so will drop these
# species_to_keep <- species_to_keep[!species_to_keep %in% c("Megachile comata", "Megachile subexilis")]
# 
# phylo_tips <- phylo$tip.label
# #only keep tips that match our species
# phylo <- ape::keep.tip(phylo, species_to_keep[species_to_keep %in% phylo_tips])
# 
phylo_matrix <- ape::vcv.phylo(phylo)
# 
# ## dropping species not in the phylogeny from the dataset
# drop.species <- unique(spec.net$GenusSpecies[!(spec.net$GenusSpecies %in% rownames(phylo_matrix))])
# 
# # here is where weights are all 0 TODO: fix
# spec.net <- spec.net[!spec.net$GenusSpecies %in% drop.species,]


## check which individuals don't have microbe data
drop.PD.NA <- unique(spec.net$UniqueID[spec.net$WeightsMicrobe == 1 &
                                         is.na(spec.net$PD)])

## drop individuals that had parasite screen but not microbe
## filling in zeros for NAs in PD
spec.net <- spec.net[!(spec.net$UniqueID %in% drop.PD.NA),] %>%
  mutate(PD = ifelse(!is.na(PD), PD, 0))



##adding abundance weights column
abund_csv <- data.frame(read.csv("../../data/sp_year_site_round.csv"))

#join abundance csv
spec.net <- join(spec.net, abund_csv)

#genus.microbes <- spec.microbes[spec.microbes$Genus == this_genus, ]

## TODO make this a function?

## what I want:
## for subset -- should be ones and zeroes ones for bombus that had microbe screening
## for weights -- log abundance weights for obligate and transient, but only for that genus 

## make genus specific weights 0s and  bombus
spec.net$BombusWeights = ifelse(spec.net$Genus=='Bombus'&spec.net$WeightsMicrobe==1, 1, 0)

## make genus specific weights 0s and 1s apis
spec.net$ApisWeights = ifelse(spec.net$Genus=='Apis'&spec.net$WeightsMicrobe==1, 1, 0)

## make genus specific weights 0s and 1s melissodes
spec.net$MelissodesWeights = ifelse(spec.net$Genus=='Melissodes'&spec.net$WeightsMicrobe==1, 1, 0)

## multiply weightspar by abundance to get abundance weights
spec.net$WeightsAbund <- spec.net$WeightsMicrobe * spec.net$AbundanceSYR

## log transformed weights abund
spec.net$LogWeightsAbund <- log(spec.net$WeightsAbund + 1)

## log weights by abundance for obligate microbes only
spec.net$LogWeightsObligateAbund <- spec.net$WeightsObligateMicrobe * spec.net$LogWeightsAbund

## log weights by abundance for transient microbes only
spec.net$LogWeightsTransientAbund <- spec.net$WeightsTransientMicrobe * spec.net$LogWeightsAbund

## log weights by abundance for obligate Bombus microbes only
spec.net$BombusLogWeightsObligateAbund <- ifelse(spec.net$Genus=='Bombus'&spec.net$LogWeightsObligateAbund>0, spec.net$LogWeightsObligateAbund, 0)

## log weights by abundance for transient Bombus microbes only
spec.net$BombusLogWeightsTransientAbund <- ifelse(spec.net$Genus=='Bombus'&spec.net$LogWeightsTransientAbund>0, spec.net$LogWeightsTransientAbund, 0)

## log weights by abundance for obligate Apis microbes only
spec.net$ApisLogWeightsObligateAbund <- ifelse(spec.net$Genus=='Apis'&spec.net$LogWeightsObligateAbund>0, spec.net$LogWeightsObligateAbund, 0)

## log weights by abundance for transient Apis microbes only
spec.net$ApisLogWeightsTransientAbund <- ifelse(spec.net$Genus=='Apis'&spec.net$LogWeightsTransientAbund>0, spec.net$LogWeightsTransientAbund, 0)

## log weights by abundance for obligate Melissodes microbes only
spec.net$MelissodesLogWeightsObligateAbund <- ifelse(spec.net$Genus=='Melissodes'&spec.net$LogWeightsObligateAbund>0, spec.net$LogWeightsObligateAbund, 0)

## log weights by abundance for transient Melissodes microbes only
spec.net$MelissodesLogWeightsTransientAbund <- ifelse(spec.net$Genus=='Melissodes'&spec.net$LogWeightsTransientAbund>0, spec.net$LogWeightsTransientAbund, 0)

## TODO add to init if works
spec.net$WeightsObligateBombus = spec.net$WeightsObligateMicrobe*spec.net$BombusWeights

## TODO add to init if works
spec.net$WeightsTransientBombus = spec.net$WeightsTransientMicrobe*spec.net$BombusWeights

## TODO add to init if works

# use pd obligate with skew normal
spec.net$WeightsObligateApis = spec.net$WeightsObligateMicrobe*spec.net$ApisWeights

#use pd transient with skew normal?
spec.net$WeightsTransientApis = spec.net$WeightsTransientMicrobe*spec.net$ApisWeights

# use pd obligate lof with skew normal?
spec.net$WeightsObligateMelissodes = spec.net$WeightsObligateMicrobe*spec.net$MelissodesWeights
# use pd transient log with gaussian or student t?
spec.net$WeightsTransientMelissodes = spec.net$WeightsTransientMicrobe*spec.net$MelissodesWeights




## 8/21/24 I think this can be done on all the spec.net data since I set up the weights differently above 

spec.net$PD.obligate.log <- log(spec.net$PD.obligate + 1)
spec.net$PD.obligate.log <- ifelse(is.na(spec.net$PD.obligate.log), 0, spec.net$PD.obligate.log)

spec.net$PD.transient.log <- log(spec.net$PD.transient + 1)
spec.net$PD.transient.log <- ifelse(is.na(spec.net$PD.transient.log), 0, spec.net$PD.transient.log)


save(spec.net, file="../../data/spec_microbes.Rdata")
