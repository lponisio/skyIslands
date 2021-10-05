## goal = input NCBI and RDP calls, make some prelim calls to help manual calling go faster
## end up with a final dataset called taxonomyRBCL.txt which we will turn into a .qza


## Use these rules
## 1	If genus matches on NBCI and RDP, use that genus.
## 2  If species matches on NCBI and RDP, use that species
## 3	If genus/sp do not match, look up samples in veg data. If one sample is present in the veg data and other is absent, use that genus
## 4	if both are present in veg data, use NCBI unless the RDP classifier is over .75 confidence, in which case, use your best judgement
## 5	If neither are present in veg, look at at top 10 NCBI BLAST hits. If any hit is in the veg data or in the region, use it as final call
## 6  If top 10 hits don't get you closer, use family level ID if it's the same for RDP/NCBI and feasible for the region
## 7 	 We have data for what bees were caught on. There are a couple species that bees were caught on, but their pollen rbcL matched with something else as the top hit. 
##       What they had been caught on did match in NCBI though, with only a single snp difference causing the known species to not be the first hit (e.g., 176/179 bp matches for Ipomoea wrightii vs 175/179 for Convolvulus arvensis). 
##       Following this logic, we manually changed I. wrightii to Con arv, and Arctotheca calendula to H. annuus.   


rm(list=ls())
ncbi <- read.table("~/Dropbox/sunflower_saved/ffar_pipeline/merged/RBCL/rbcl_classified_NCBI.txt",
                   sep="\t", header=TRUE)


rdp <- read.table("~/Dropbox/sunflower_saved/ffar_pipeline/merged/RBCL/rbcl_classified_rdp.txt",
                  sep="\t", header=TRUE)

## rdp cleaning
rdp$clean_family <- gsub("f__", "", rdp$family)
rdp$clean_family <- sapply(strsplit(rdp$clean_family, "_"),
                           function(x) x[1])

rdp$clean_genus <- gsub("g__", "", rdp$genus)
rdp$clean_genus <- sapply(strsplit(rdp$clean_genus, "_"),
                          function(x) x[1])

rdp$clean_species <- gsub("s__", "", rdp$species)
rdp$GenusSpecies <- sapply(strsplit(rdp$clean_species, "_"),
                           function(x) x[1])

## ncbi cleaning

ncbi$genus <- sapply(strsplit(ncbi$name, "[ ]"),
                     function(x) x[1])

ncbi$species <- sapply(strsplit(ncbi$name, "[ ]"),
                       function(x) x[2])
ncbi$GenusSpecies <- paste(ncbi$genus, ncbi$species)


## plants IDed by hand at sites
plants <- read.csv("~/Dropbox/sunflower_saved/data/raw/plants.csv")

veg.genera <- unique(plants$Genus)
veg.genera <- veg.genera[!is.na(veg.genera)]

veg.sp <- unique(plants$Species)
veg.sp <- veg.sp[!is.na(veg.sp)]

nrow(ncbi)
nrow(rdp)

ids <- data.frame(Sample=rdp$ID,
                  GenusSpeciesRDP=rdp$GenusSpecies,
                  GenusSpeciesRDPMatch=rdp$species_match,
                  GenusRDP=rdp$clean_genus,
                  GenusRDPMatch=rdp$genus_match)	       

ids$GenusSpeciesNCBI <- ncbi$GenusSpecies[match(ids$Sample,
                                                ncbi$Query)]
ids$GenusNCBI <- ncbi$genus[match(ids$Sample,
                                  ncbi$Query)]

ids$GenusRDP_NCBI_match <- ids$GenusNCBI == ids$GenusRDP
ids$SpeciesRDP_NCBI_match <- ids$GenusSpeciesNCBI == ids$GenusSpeciesRDP

ids$FinalGenusSpecies <- NA
ids$FinalGenus <- NA


## genera that match between RDP and NCBI
ids$FinalGenus[ids$GenusRDP_NCBI_match & !is.na(ids$GenusRDP_NCBI_match)]  <- ids$GenusRDP[ids$GenusRDP_NCBI_match & !is.na(ids$GenusRDP_NCBI_match)]

## and are also in the veg data IDed by hand
unique(ids$FinalGenus)
ids$FinalGenusInVeg <- NA
ids$FinalGenusInVeg[ids$FinalGenus %in% veg.genera] <- 1
ids$FinalGenusInVeg[!ids$FinalGenus %in% veg.genera & !is.na(ids$FinalGenus)] <- 0
unique(ids$FinalGenus[ids$FinalGenusInVeg == 0])

## make a final call 

## IDs that have a solid match but are not in the hand collected veg dats 
unique(ids$FinalGenus[ids$FinalGenusInVeg == 0][ids$GenusRDPMatch > 0.8])
## 

## IDs that have an okay match but are not in the hand collected veg dats 
unique(ids$FinalGenus[ids$FinalGenusInVeg == 0][ids$GenusRDPMatch >0.6 & ids$GenusRDPMatch < 0.8])


acceptable <- c("Xxxx")

not.possible <- c("Xxxx")

## if acceptable, gets final genus check
ids$FinalGenusConfirmed  <- NA

## confirmed if matches veg
ids$FinalGenusConfirmed[ids$FinalGenusInVeg == 1] <- TRUE

## confirmed if hand checked
ids$FinalGenusConfirmed[ids$FinalGenus %in% acceptable] <- TRUE

ids$FinalGenus[ids$FinalGenus %in% not.possible] <- paste("Unknown near", ids$FinalGenus[ids$FinalGenus %in% not.possible])

ids$FinalGenusConfirmed[ids$FinalGenus %in% not.possible] <- TRUE

## ************************************************************************************
## species that match between RDP and NCBI
ids$FinalGenusSpecies[ids$SpeciesRDP_NCBI_match & !is.na(ids$SpeciesRDP_NCBI_match)]  <- ids$GenusSpeciesRDP[ids$SpeciesRDP_NCBI_match & !is.na(ids$SpeciesRDP_NCBI_match)]

## and are also in the veg data IDed by hand
unique(ids$FinalGenusSpecies)
ids$FinalGenusSpeciesInVeg <- NA
ids$FinalGenusSpeciesInVeg[ids$FinalGenusSpecies %in% veg.sp] <- 1
ids$FinalGenusSpeciesInVeg[!ids$FinalGenusSpecies %in% veg.sp & !is.na(ids$FinalGenusSpecies)] <- 0

unique(ids$FinalGenusSpecies[ids$FinalGenusSpeciesInVeg == 0])
unique(ids$FinalGenusSpecies[ids$FinalGenusSpeciesInVeg == 0][ids$GenusRDPMatch >0.8])

## "Vicia cirrhosa"        "Chlorophytum nimmonii"  "Rubus idaeus"
## none of these are at all likely

## where NCBI didn't have a match, if RDP is > 50%
ids$FinalGenus[is.na(ids$GenusRDP_NCBI_match) & ids$GenusRDPMatch > 0.5] <- ids$GenusRDP[is.na(ids$GenusRDP_NCBI_match) & ids$GenusRDPMatch > 0.5]
## this never happens

not.possible <- c("Xxxx")
####

## if acceptable, gets final species check
ids$FinalGenusSpeciesConfirmed  <- NA

## confirmed if matches veg
ids$FinalGenusSpeciesConfirmed[ids$FinalGenusSpeciesInVeg == 1] <- TRUE


## confirmed if hand checked

ids$FinalGenusSpecies[ids$FinalGenusSpecies %in% not.possible] <- paste("Unknown near", ids$FinalGenusSpecies[ids$FinalGenusSpecies %in% not.possible])

ids$FinalGenusSpeciesConfirmed[ids$FinalGenusSpecies %in% not.possible] <- TRUE

ids$FinalGenusSpeciesInVeg <- NULL
ids$FinalGenusInVeg <- NULL

write.csv(ids, "~/Dropbox/sunflower_saved/ffar_pipeline/merged/RBCL/taxonomyRBCLworkbook.csv", row.names=FALSE)


## editing final csv
final <- read.csv("~/Dropbox/sunflower_saved/ffar_pipeline/merged/RBCL/taxonomyRBCLworkbook.csv",
                  sep="\t", header=TRUE)


