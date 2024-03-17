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
library(scales)
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

## FULL MICROBE DATASET

## pull out 16s columns
microbes <- colnames(spec.net)[grepl("16s:", colnames(spec.net))] 
microbes <- microbes[microbes %in% tree.16s$tip.label]

## check which only have NAs in these columns (not screened) and drop them
screened.microbes <- apply(spec.net, 1, function(x) all(is.na(x[microbes])))
spec.microbes <- spec.net[!screened.microbes, ]

#3 rows have 0 for all microbes, need to drop
spec.microbes <- spec.microbes[rowSums(spec.microbes[,microbes])!=0,]


## QUESTION: should include root = TRUE? if false gives warning 3x
## warning: Rooted tree and include.root=TRUE argument required to calculate PD of single-species communities. Single species community assigned PD value of NA.

## phylogenetic distance function, modified from picante
PD <- apply(spec.microbes[,microbes], 1, function(x){
  this.bee <- x[x > 0]
  this.tree <- prune.sample(t(this.bee), tree.16s)
  picante::pd(t(this.bee), this.tree, include.root = TRUE)
  #browser()
})

PD <- do.call(rbind, PD)
spec.microbes <- cbind(spec.microbes, PD)

## Merge back onto specimen data
spec.net <- merge(spec.net, spec.microbes, all.x=TRUE)

## change microbe NAs to 0
spec.net <- spec.net %>%
  mutate_at(vars(all_of(microbes)), ~replace_na(PD, 0))
  
spec.net[,microbes][is.na(spec.net[,microbes])] <- 0


## ONLY OBLIGATE MICROBES DATASET


## splitting out obligate bee microbes based on Zheng and Moran paper
bee.obligates <- "Lactobacillus|Bifidobacterium|Snodgrassella|Gilliamella|Frischella|Bartonella|Commensalibacter"

## this is a list of the microbe strains that contain the known obligate bee microbe genera
bee.obligate.microbes <- microbes[grepl(bee.obligates, microbes, fixed=FALSE)]
bee.obligate.microbes <- bee.obligate.microbes[bee.obligate.microbes %in% tree.16s$tip.label]
## now need to subset spec.microbes to be just the microbe columns with bee obligates and calculate PD

## phylogenetic distance function, modified from picante
PD.obligate <- apply(spec.microbes[,bee.obligate.microbes], 1, function(x){
  tryCatch({
    this.bee <- x[x > 0]
    this.tree <- prune.sample(t(this.bee), tree.16s)
    pd_value <- picante::pd(t(this.bee), this.tree, include.root = TRUE)
    if (is.null(pd_value) || length(pd_value) == 0) pd_value <- 0  # Assign zero if PD is NULL or empty list
    else pd_value[[1]]  # Extract the first element if PD is a list
    #browser()
    # Return PD and SR as a list
    #data.frame(PD = pd_value[[1]], SR = pd_value[[2]])
  }, error = function(e) {
    if (grepl("Tree has no branch lengths", e$message)) {
      # If error is due to tree having no branch lengths, return zero for PD and SR
      return(list(PD = 0, SR = 0))
    } else {
      # If it's a different error, re-raise it
      stop(e)
    }
  })
})

# Convert the result into a dataframe
result_df <- do.call(rbind, PD.obligate)
rownames(result_df) <- NULL  # Remove row names

# Optionally, you can rename the columns
colnames(result_df) <- c("PD", "SR")

PD <- do.call(rbind, PD)
spec.microbes <- cbind(spec.microbes, PD)

## Merge back onto specimen data
spec.net <- merge(spec.net, spec.microbes, all.x=TRUE)

## change microbe NAs to 0
spec.net <- spec.net %>%
  mutate_at(vars(all_of(microbes)), ~replace_na(PD, 0))

spec.net[,microbes][is.na(spec.net[,microbes])] <- 0



## create a dumby varaible "Weight" to deal with the data sets being at
## different levels to get around the issue of having to pass in one
## data set into brms
spec.net$YearSR <- paste(spec.net$Year, spec.net$SampleRound, sep=";")
spec.net$YearSRGenusSpecies <- paste(spec.net$YearSR, spec.net$GenusSpecies, sep=";")

## will need to modify when we have multiple years
spec.net <- makeDataMultiLevel(spec.net, "Site", "YearSR")

#the 6 species getting dropped are still present at this step

if(make.plots == FALSE){
  spec.net[, variables.to.log] <- log(spec.net[,variables.to.log])
}
#the 6 species getting dropped are still present at this step

##  center all of the x variables, need to use unique values to avoid
##  repetition by the number of specimens

if(make.plots == FALSE){
spec.net <- standardizeVars(spec.net, vars_yearsr, "YearSR")

#the 6 species getting dropped are still present at this step

spec.net <- standardizeVars(spec.net, vars_sp, "YearSRGenusSpecies")
}
##load tree from :
##Henriquez Piskulich, Patricia Andrea; Hugall, Andrew F.; Stuart-Fox, Devi (2023). A supermatrix phylogeny of the world’s bees (Hymenoptera: Anthophila) [Dataset]. Dryad. https://doi.org/10.5061/dryad.80gb5mkw1

phylo <- ape::read.tree("../data/BEE_mat7_fulltree.nwk")


## i think we need to drop the tips that aren't in our dataset
## changing sp labels to match phylogeny
# spec.net$GenusSpecies <- gsub("Megachile gemula fulvogemula", "Megachile gemula", spec.net$GenusSpecies)
# spec.net$GenusSpecies <- gsub("Megachile melanophaea rohweri", "Megachile melanophaea", spec.net$GenusSpecies)


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

#the 6 species getting dropped are still present at this step

phylo_tips <- phylo$tip.label
#only keep tips that match our species
phylo <- ape::keep.tip(phylo, species_to_keep[species_to_keep %in% phylo_tips])

phylo_matrix <- ape::vcv.phylo(phylo)

## dropping species not in the phylogeny from the dataset
drop.species <- unique(spec.net$GenusSpecies[!(spec.net$GenusSpecies %in% rownames(phylo_matrix))])
spec.net <- spec.net[!spec.net$GenusSpecies %in% drop.species,]

## create a dumby varaible "WeightPar" for the parasite data. The
## original intention was to keep stan from dropping data for
## site-level models, but weight is 0 for parasite models.
spec.net <- prepParasiteWeights()

#the 6 species getting dropped are still present at this step

#genus_pd_fit <- function(spec.net, this_genus, num_iter){



#the 6 species getting dropped are still present at this step

## check which individuals don't have microbe data
drop.PD.NA <- unique(spec.net$UniqueID[spec.net$WeightsMicrobe == 1 &
                                         is.na(spec.net$PD)])

## drop individuals that had parasite screen but not microbe
## filling in zeros for NAs in PD
spec.net <- spec.net[!(spec.net$UniqueID %in% drop.PD.NA),] %>%
  mutate(PD = ifelse(!is.na(PD), PD, 0))

#the 6 species getting dropped are still present at this step

##adding abundance weights column
abund_csv <- data.frame(read.csv("../../data/sp_year_site_round.csv"))

#join abundance csv
spec.net <- join(spec.net, abund_csv)


#the 6 species getting dropped are still present at this step

#genus.microbes <- spec.microbes[spec.microbes$Genus == this_genus, ]



#multiply weightspar by abundance to get abundance weights
spec.net$WeightsAbund <- spec.net$WeightsMicrobe * spec.net$AbundanceSYR
spec.net$LogWeightsAbund <- log(spec.net$WeightsAbund + 1)
spec.net$BombusWeights <- ifelse(spec.net$Genus=='Bombus'&spec.net$LogWeightsAbund>0, spec.net$LogWeightsAbund, 0)
spec.net$ApisWeights <- ifelse(spec.net$Genus=='Apis'&spec.net$LogWeightsAbund>0, spec.net$LogWeightsAbund, 0)
spec.net$MelissodesWeights <- ifelse(spec.net$Genus=='Melissodes'&spec.net$LogWeightsAbund>0, spec.net$LogWeightsAbund, 0)



spec.all <- spec.net

## bombus only data
spec.bombus <- spec.all
#multiply weightsMicrobe by abundance to get abundance weights
spec.bombus$WeightsAbund <- spec.bombus$WeightsMicrobe * spec.bombus$AbundanceSYR
spec.bombus$LogWeightsAbund <- log(spec.bombus$WeightsAbund + 1)
spec.bombus$LogWeightsAbund[spec.bombus$Genus != "Bombus"] <- 0
spec.bombus$WeightsAbund[spec.bombus$Genus != "Bombus"] <- 0

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
#multiply weightsMicrobe by abundance to get abundance weights
spec.apis$WeightsAbund <- spec.apis$WeightsMicrobe * spec.apis$AbundanceSYR
spec.apis$LogWeightsAbund <- log(spec.apis$WeightsAbund + 1)
spec.apis$LogWeightsAbund[spec.apis$Genus != "Apis"] <- 0
spec.apis$WeightsAbund[spec.apis$Genus != "Apis"] <- 0

## melissodes only data
spec.melissodes <- spec.all
#multiply weightsMicrobe by abundance to get abundance weights
spec.melissodes$WeightsAbund <- spec.melissodes$WeightsMicrobe * spec.melissodes$AbundanceSYR
spec.melissodes$LogWeightsAbund <- log(spec.melissodes$WeightsAbund + 1)
spec.melissodes$LogWeightsAbund[spec.melissodes$Genus != "Melissodes"] <- 0
spec.melissodes$WeightsAbund[spec.melissodes$Genus != "Melissodes"] <- 0

# spec.apidae <- spec.all
# spec.apidae$WeightsMicrobe[spec.melissodes$Family != "Apidae"] <- 0

### troubleshooting
these_are_getting_dropped <- c("Anthidium mormonum","Bombus centralis","Bombus huntii","Dufourea maura","Melissodes confusus")
#View(spec.net[spec.net$GenusSpecies %in% these_are_getting_dropped,])
