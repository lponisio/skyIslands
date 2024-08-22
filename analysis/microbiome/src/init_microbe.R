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


########## 

## weights by abundance for transient microbes only
# spec.net$WeightsTransientAbund <- spec.net$WeightsTransientMicrobe * spec.net$AbundanceSYR
# 
# ## make weights for each bombus separately for obligate and transient
# spec.net$BombusObligateWeights <- ifelse(spec.net$Genus=='Bombus'&spec.net$WeightsObligateAbund>0, spec.net$WeightsObligateAbund, 0)
# spec.net$LogBombusObligateWeightsAbund <- log(spec.net$WeightsObligateAbund + 1)
# spec.net$BombusTransientWeights <- ifelse(spec.net$Genus=='Bombus'&spec.net$WeightsTransientAbund>0, spec.net$WeightsTransientAbund, 0)
# 
# ## make weights for each genus separately for obligate and transient
# spec.net$ApisObligateWeights <- ifelse(spec.net$Genus=='Apis'&spec.net$WeightsObligateAbund>0, spec.net$WeightsObligateAbund, 0)
# 
# ## make weights for each genus separately for obligate and transient
# spec.net$ApisTransientWeights <- ifelse(spec.net$Genus=='Apis'&spec.net$WeightsTransientAbund>0, spec.net$WeightsTransientAbund, 0)
# 
# ## make weights for each genus separately for obligate and transient
# spec.net$MelissodesObligateWeights <- ifelse(spec.net$Genus=='Melissodes'&spec.net$WeightsObligateAbund>0, spec.net$WeightsObligateAbund, 0)

## make weights for each genus separately for obligate and transient
# spec.net$MelissodesTransientWeights <- ifelse(spec.net$Genus=='Melissodes'&spec.net$WeightsTransientAbund>0, spec.net$WeightsTransientAbund, 0)

## TODO check if this can be removed now

# #multiply weightsMicrobe by abundance to get abundance weights
# spec.bombus$WeightsAbund <- spec.bombus$WeightsMicrobe * spec.bombus$AbundanceSYR
# spec.bombus$WeightsObligateAbund <- spec.bombus$WeightsObligateMicrobe * spec.bombus$AbundanceSYR
# spec.bombus$LogWeightsObligateAbund <- log(spec.bombus$WeightsObligateAbund + 1)
# spec.bombus$LogWeightsObligateAbund[spec.bombus$Genus != "Bombus"] <- 0



# spec.net$ApisWeights <- ifelse(spec.net$Genus=='Apis'&spec.net$LogWeightsAbund>0, spec.net$LogWeightsAbund, 0)
# spec.net$MelissodesWeights <- ifelse(spec.net$Genus=='Melissodes'&spec.net$LogWeightsAbund>0, spec.net$LogWeightsAbund, 0)


# 
# spec.all <- spec.net
# 
# ## bombus only data
# spec.bombus <- spec.all
# #multiply weightsMicrobe by abundance to get abundance weights
# spec.bombus$WeightsAbund <- spec.bombus$WeightsMicrobe * spec.bombus$AbundanceSYR
# spec.bombus$WeightsObligateAbund <- spec.bombus$WeightsObligateMicrobe * spec.bombus$AbundanceSYR
# spec.bombus$LogWeightsObligateAbund <- log(spec.bombus$WeightsObligateAbund + 1)
# spec.bombus$LogWeightsObligateAbund[spec.bombus$Genus != "Bombus"] <- 0
# spec.bombus$LogWeightsAbund <- log(spec.bombus$WeightsAbund + 1)
# spec.bombus$LogWeightsAbund[spec.bombus$Genus != "Bombus"] <- 0
# spec.bombus$WeightsAbund[spec.bombus$Genus != "Bombus"] <- 0
# spec.bombus$WeightsTransientAbund <- spec.bombus$WeightsTransientMicrobe * spec.bombus$AbundanceSYR
# spec.bombus$LogWeightsTransientAbund <- log(spec.bombus$WeightsTransientAbund + 1)
# spec.bombus$LogWeightsTransientAbund[spec.bombus$Genus != "Bombus"] <- 0
# 
# # 
# # ## megachile comate and megachile subexilis are not in phylogeny so will drop these
# # species_to_keep <- species_to_keep[!species_to_keep %in% c("Megachile comata", "Megachile subexilis")]
# # 
# # phylo_tips <- phylo$tip.label
# # #only keep tips that match our species
# # phylo <- ape::keep.tip(phylo, species_to_keep[species_to_keep %in% phylo_tips])
# # 
# # phylo_matrix <- ape::vcv.phylo(phylo)
# 
# 
# ## apis only data
# spec.apis <- spec.all
# #multiply weightsMicrobe by abundance to get abundance weights
# spec.apis$WeightsAbund <- spec.apis$WeightsMicrobe * spec.apis$AbundanceSYR
# spec.apis$LogWeightsAbund <- log(spec.apis$WeightsAbund + 1)
# spec.apis$LogWeightsAbund[spec.apis$Genus != "Apis"] <- 0
# spec.apis$WeightsAbund[spec.apis$Genus != "Apis"] <- 0
# spec.apis$WeightsObligateAbund <- spec.apis$WeightsObligateMicrobe * spec.apis$AbundanceSYR
# spec.apis$LogWeightsObligateAbund <- log(spec.apis$WeightsObligateAbund + 1)
# spec.apis$LogWeightsObligateAbund[spec.apis$Genus != "Apis"] <- 0
# spec.apis$WeightsTransientAbund <- spec.apis$WeightsTransientMicrobe * spec.apis$AbundanceSYR
# spec.apis$LogWeightsTransientAbund <- log(spec.apis$WeightsTransientAbund + 1)
# spec.apis$LogWeightsTransientAbund[spec.apis$Genus != "Apis"] <- 0
# 
# ## melissodes only data
# spec.melissodes <- spec.all
# 
# #multiply weightsMicrobe by abundance to get abundance weights
# spec.melissodes$WeightsAbund <- spec.melissodes$WeightsMicrobe * spec.melissodes$AbundanceSYR
# spec.melissodes$LogWeightsAbund <- log(spec.melissodes$WeightsAbund + 1)
# spec.melissodes$LogWeightsAbund[spec.melissodes$Genus != "Melissodes"] <- 0
# spec.melissodes$WeightsAbund[spec.melissodes$Genus != "Melissodes"] <- 0
# spec.melissodes$WeightsObligateAbund <- spec.melissodes$WeightsObligateMicrobe * spec.melissodes$AbundanceSYR
# spec.melissodes$LogWeightsObligateAbund <- log(spec.melissodes$WeightsObligateAbund + 1)
# spec.melissodes$LogWeightsObligateAbund[spec.melissodes$Genus != "Melissodes"] <- 0
# spec.melissodes$WeightsTransientAbund <- spec.melissodes$WeightsTransientMicrobe * spec.melissodes$AbundanceSYR
# spec.melissodes$LogWeightsTransientAbund <- log(spec.melissodes$WeightsTransientAbund + 1)
# spec.melissodes$LogWeightsTransientAbund[spec.melissodes$Genus != "Melissodes"] <- 0

# spec.apidae <- spec.all
# spec.apidae$WeightsMicrobe[spec.melissodes$Family != "Apidae"] <- 0

## 6-3-24 RH removing the +1 since we are dropping zeros anyway
## 6-4-24 RH leaving +1 in because logging the tiny numbers is making the distribution even weirder
## which seemingly negates the point of transforming the data... 


## 8/21/24 I think this can be done on all the spec.net data since I set up the weights differently above 

spec.net$PD.obligate.log <- log(spec.net$PD.obligate + 1)
spec.net$PD.obligate.log <- ifelse(is.na(spec.net$PD.obligate.log), 0, spec.net$PD.obligate.log)

spec.net$PD.transient.log <- log(spec.net$PD.transient + 1)
spec.net$PD.transient.log <- ifelse(is.na(spec.net$PD.transient.log), 0, spec.net$PD.transient.log)

# spec.apis$PD.obligate.log <- log(spec.apis$PD.obligate + 1)
# spec.apis$PD.obligate.log <- ifelse(is.na(spec.apis$PD.obligate.log), 0, spec.apis$PD.obligate.log)
# 
# spec.apis$PD.transient.log <- log(spec.apis$PD.transient + 1)
# spec.apis$PD.transient.log <- ifelse(is.na(spec.apis$PD.transient.log), 0, spec.apis$PD.transient.log)
# 
# spec.melissodes$PD.obligate.log <- log(spec.melissodes$PD.obligate + 1)
# spec.melissodes$PD.obligate.log <- ifelse(is.na(spec.melissodes$PD.obligate.log), 0, spec.melissodes$PD.obligate.log)
# 
# spec.melissodes$PD.transient.log <- log(spec.melissodes$PD.transient + 1)
# spec.melissodes$PD.transient.log <- ifelse(is.na(spec.melissodes$PD.transient.log), 0, spec.melissodes$PD.transient.log)


save(spec.net, file="../../data/spec_microbes.Rdata")
