## *******************************************************
## Pseudocode based on Trojelsgaard paper - KCT 05/2018
## *******************************************************
## setwd("~/Dropbox/skyIslands")
rm(list=ls())
setwd("~/Dropbox/skyIslands/analysis")
library(gtools)
library(vegan)
load('../data/spec.Rdata')
source("../dataPrep/src/misc.R")

#Prepare networks
prep.dis.mat.pol <- samp2site.spp(site=spec[,"Site"],
                              spp=spec[,"GenusSpecies"],
                              abund=rep(1, nrow(spec)),
                              FUN=sum)


prep.dis.mat.plant <- samp2site.spp(site=spec[,"Site"],
                              spp=spec[,"PlantGenusSpecies"],
                              abund=rep(1, nrow(spec)),
                              FUN=sum)

prep.dis.mat.int <- samp2site.spp(site=spec[,"Site"],
                                  spp=spec[,"Int"],
                                  abund=rep(1,nrow(spec)),
                                  FUN=sum)

#Run dissimilarity index vegdist()
dis.mat.pol <- vegdist(prep.dis.mat.pol, "jaccard")
dis.mat.plant <- vegdist(prep.dis.mat.plant, "jaccard")
dis.mat.int <- vegdist(prep.dis.mat.int, "jaccard")

#Run null model based on 1000 random assignments of species and interactions to networks according to a probability distribution derived from their actual occurrences across the sites (widespread species more likely to be drawn; numbers of species and interactions assigned to a random network constrained to equal empirical numbers)
		##Calculate Sorenson similarity (above formula) for each random assignment - calculate mean beta and SD
	##Calculate deviation from randomness using z-scores: (beta(empirical) - mean(beta(resampled)))/SD(resampled)
		##mean(beta(resampled)) and SD(resampled) = mean and SD of the similarities achieved from the 1000 random assortments of species and interactions
			##positive z-score suggests two communities are more similar than random; vice versa for negative z-score
	##Compare empirical Sorensen similarity values and derived z-scores with geographic distance using Mantel tests with 1000 permutations
		##VEGAN v. 2.0-8 package for R


##2. Explain the pattern of dissimilarity among networks by analyzing the underlying interaction turnover.
	##Estimate the relative contribution of species-driven interaction turnover and rewiring
		
    ##I(rewired) = number of interactions that change between shared species of A and B (A and B are two networks)
#Sounds like an %in% situation to me!
SI.agg <- aggregate(spec$GenusSpecies, list(site=spec$Site, plant=spec$PlantGenusSpecies, pol=spec$GenusSpecies), length)
SI.split<-split(SI.agg, SI.agg[[1]], drop=TRUE)

Int.agg <- aggregate(spec$Int, list(site=spec$Site, int=spec$Int), length)
Int.split <- split(Int.agg, Int.agg[[1]], drop=TRUE)

list.of.plants <- list()
for(i in 1:length(SI.split)){
  this.site <- SI.split[[i]]
  this.site.plants <- unique(this.site$plant)
  list.of.plants[[i]] <- this.site.plants
}
names(list.of.plants) <- names(SI.split)

list.of.pols <- list()
for(i in 1:length(SI.split)){
  this.site <- SI.split[[i]]
  this.site.pols <- unique(this.site$pol)
  list.of.pols[[i]] <- this.site.pols
}
names(list.of.pols) <- names(SI.split)

list.of.ints <- list()
for(i in 1:length(Int.split)){
  this.site <- Int.split[[i]]
  this.site.ints <- unique(this.site$int)
  list.of.ints[[i]] <- this.site.ints
}
names(list.of.ints) <- names(Int.split)

##For some reason we haven't eliminated all blank plants....need to figure that out... - KCT 6/11/18

		##I(species) = number of interactions that change due to changes in species composition
			##I(rewired) + I(species) = total number of different interactions between A and B
			##proportion of turnover due to rewiring: I(rewired)/(I(rewired) + I(species))
			##proportion of turnover due to species turnover: I(species)/(I(rewired)+I(species))
	##Relate the above two trends with geographical distance
	##Look at drivers of species-driven interaction turnover (I(species))
		##I(pla) = number of interactions between non-shared plants and shared pollinators
		##I(pol) = number of interactions between non-shared pollinators and shared plant species
		##I(pol+pla) = number of interactions between non-shared plants and non-shared pollinators (complete turnover)
			##I(species) = I(pol) + I(pla) + I(pol+pla)
			##Fraction of I(species) explained by replacement of pollinators: I(pol)/I(species)
			##Fraction of I(species) explained by replacement of plants: I(pla)/I(species)
			##Fraction of I(species) explained by replacement of both: I(pol+pla)/I(species)
	##Something about cubic smoothed splines...not quite sure what this is...


##3. Quantify partner fidelity by comparing actual partner composition of individual plant and pollinator species with a random set of partners - explore ecological and geographical determinants of partner fidelity
	##Calculate average Sorenson similarity (Sor(emp)) in partner composition for all plants and pollinators occurring in more than two networks
	##Randomly select same number of partners for that species from all the available species at the site and recalculate average Sorenson similarity between all pairwise comparisons (repeat 1000 times, calculate z-scores as (Sor(emp) - Sor(rand))/Sor(SD))
		##Resolved partners of plant species (pollinators) to species, genus, family, and order
		##Resolved partners of pollinator species (plants) to species and family

	##Assess how partner similarities were affected by geographic distance between networks, difference in abundance of focal species between compared networks, difference in number of interaction partners of focal species between compared networks, number of sites at which focal species existed
		##Linear mixed-effects model and model averaging
			##Species-specific partner similarity values as a response variable
			##Species as a random factor
			##Random intercepts and random slopes allowed
			##Run models with all possible combinations of fixed parameters (geographic distance, change in abundance, change in linkage levele, geographic distribution)
				##Perform a model averaging procedure of the parameter estimates using the Akaike weights (not sure what this means...)
			##Use overall similarity of the partner community as a fixed explanatory variable
				##LME4 v. 1.0-5 and MUMIN v. 1.9.13 for linear mixed models and model averaging, respectively
