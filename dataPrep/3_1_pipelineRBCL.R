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
setwd("~/")
source("lab_paths.R")
local.path
setwd(local.path)
setwd("skyIslands_saved")
library(dplyr)

ncbi <- read.csv("SI_pipeline/merged/RBCL/rbcl_classified_NCBI.csv")

rdp <- read.table("SI_pipeline/merged/RBCL/rbcl_classified_rdp.txt",
                  sep="\t", header=FALSE)

## loading old rdp run to pull column names from
rdp_old <- read.table("SI_pipeline/merged/RBCL/rbcl_classified_rdp_old.txt",
                  sep="\t", header=TRUE)

old_colnames <- colnames(rdp_old)

cols_to_drop <- c(7,10,13,16,19,22,25)

rdp <- rdp %>% select(!cols_to_drop)

colnames(rdp) <- old_colnames

## read.table is counting the first row as the header, we need
## to add in the right column names

## colnames ID kingdom	kingdom_match	phylum	phylum_match
## 

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

# Identify which samples were dropped from NCBI (looks like they blast to non-plant DNA)
library(arsenal)
ncbi$ID <- ncbi$Query

summary(comparedf(rdp, ncbi, by='ID'))

## plants IDed by hand at sites

plants <- read.csv(" ~/Dropbox/skyIslands_saved/data/relational/relational/traditional/veg-complete.csv")
                       
veg.genera <- unique(plants$PlantGenus)
veg.genera <- veg.genera[!is.na(veg.genera)]

veg.sp <- unique(paste(plants$PlantGenus, plants$PlantSpecies))
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

# based on USDA plant database - these genera are represented in southwestern US, but keep in mind that many of the spp hits we are getting are not found there
acceptable <- c("Pedicularis", "Conyza", "Achillea", "Allium", "Alnus", "Penstemon","Tiquilia",
                "Aloe","Amsinckia","Erigeron", "Arceuthobium", "Arctium", "Cichorium","Tanacetum",
                "Nasturtium", "Hypericum", "Campanula","Linum","Potentilla", "Cruciata", "Picea",
                "Eremogone", "Lupinus","Medicago","Oxalis","Erodium","Holodiscus","Plantago",
                "Rubus","Centaurium", "Thalictrum", "Astragalus","Bistorta","Juncus","Trifolium",
                "Prunus", "Hedyotis", "Stenaria", "Geranium","Senecio", "Brickellia","Capparis", 
                "Cercocarpus", "Amblyopappus", "Cicuta", "Zizia", "Nothoscordum", "Helianthus", 
                "Phacelia", "Veronica", "Arctotheca","Cerastium","Cichorium", "Cucurbita", "Daucus",
                "Ellisia", "Erigeron", "Ipomoea","Oxytropis", "Phacelia", "Rubus", "Scandix",
                "Sisyrinchium", "Spiranthes","Symphoricarpos", "Tamarix", "Vicia", "Ipomoea",
                "Robinia", "Galium", "Prunella", "Heuchera", "Valeriana", "Astragalus","Scorphularia",
                "Fragaria","Penstemon", "Dalea", "Calochortus", "Rosa", "Bistorta", "Swertia", "Primula",
                "Geum", "Betula", "Sedum","Garrya", "Melilotus", "Pyrola","Euphorbia", "Sanguisorba",
                "Cinnamomum","Carex","Sium","Prosopis", "Prosopis","Mentzelia","Iris","Hydrophyllum","Ziziphus",
                "Ulmus", "Chlorophytum","Eriobotrya", "Physocarpus", "Rorippa", "Polemonium", "Mimosa",
                "Crotalaria", "Glycine")

# based on USDA plant database, Discover Life, and natural history papers for some genera - these seem unlikely and should be flagged
not.possible <- c("Olearia", "Linochilus", "Acridocarpus","Sheareria","Formania", "Apodytes", 
                  "Maharanga", "Cratoxylum", "Ascolepis","Argophyllum", "Beloglottis", "Ophiopogon", #Beloglottis in FL, Hesperomannia in HI, Weigela + Liriodendron in CA
                  "Hesperomannia", "Sinocarum", "Strychnos", "Echiochilon","Indofevillea", "Enydra",
                  "Weigela", "Liriodendron")

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

## check if species is possible                       
ids$FinalGenusSpecies[ids$FinalGenusSpecies %in% not.possible] <- paste("Unknown near", ids$FinalGenusSpecies[ids$FinalGenusSpecies %in% not.possible])

ids$FinalGenusSpeciesConfirmed[ids$FinalGenusSpecies %in% not.possible] <- TRUE

ids$FinalGenusSpeciesInVeg <- NULL
ids$FinalGenusInVeg <- NULL

## merge data with sequences to do manual checks
library(dplyr)

datseq <- read.csv('~/Dropbox/skyIslands_saved/SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/Seqs_SampleNames.csv')
ID_Manual <- ids %>% left_join(datseq, by="Sample")
write.csv(ID_Manual, "~/Dropbox/skyIslands_saved/SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/ID-Manual-LJJComments.csv",row.names=FALSE) # worked on this file doing manual calls

## hand-checked samples following protocol above with (1) NCBI top 100 hits, (2) SI veg database, (3) USDA plant database, and sparingly (4) Discover Life. 
## "Conservative" calls have high confidence at that taxonomic resolution (LJJ). 
## "With uncertainty" calls have low confidence but based on veg data, USDA database, and sometimes missing seq data in NCBI for "likely" calls based on one 
                       ## of our databases (veg or USDA), one of the calls listed in that column may be correct. 
                       ## LJJ has some comments on these uncertain calls in the file "ID-Manual-LJJComments.xslx" (open with google sheets or other software that supports comments to view them).

manualcalls <- read.csv('~/Dropbox/skyIslands_saved/SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/PollenID_Manual_LJJ.csv')

## remove bacterial sequences that made it through, "bad" sequences (< 90% identity to any NCBI hits and/or no consistency at the class level for first 100 NCBI hits), 
## and any plant ID's that don't exist at those sites

manualcalls_clean <- subset(manualcalls, IDCall_Conservative != "BACTERIA" & IDCall_Conservative != "BAD" & USDA_PlantDatabase_CallNearSI != "NO")                     
                       
write.csv(manualcalls_clean, "~/Dropbox/skyIslands_saved/SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCLworkbook.csv", row.names=FALSE)

## editing final csv
final <- read.csv("~/Dropbox/skyIslands_saved/SI_pipeline/R2018/2023_sequence_results_raw/merged/RBCL/taxonomyRBCLworkbook.csv",
                  sep="\t", header=TRUE)


