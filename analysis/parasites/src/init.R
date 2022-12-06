library(ggplot2)
library(brms)
library(dplyr)


load('../../data/spec.Rdata')
site.sum <- read.csv("../../data/sitestats.csv")

parasites <- c(#"AspergillusSpp", ## problematic parasite!
               "AscosphaeraSpp",
               "ApicystisSpp",
               "CrithidiaExpoeki",
               "CrithidiaBombi",
               "NosemaBombi",
               "NosemaCeranae")

## subset to the bees we screened and the screenings worked
## spec <- spec[spec$Apidae == 1 &
##              !is.na(spec$Apidae),]

spec <- merge(spec, site.sum, all.x=TRUE)

traits <-
    read.csv("../../../skyIslands_saved/data/raw/bee_traits.csv")
traits$GenusSpecies <- fix.white.space(traits$GenusSpecies)
traits <- traits[, c("GenusSpecies", "Sociality", "Lecty", "MeanITD"),]

net.traits <- read.csv("../../data/traits.csv")
net.traits <- net.traits[, c("GenusSpecies", "r.degree"),]

traits <- merge(traits, net.traits, by="GenusSpecies", all.x=TRUE)

spec <- merge(spec, traits, all.x=TRUE, by="GenusSpecies")

dir.create(path="saved", showWarnings = FALSE)
dir.create(path="saved/tables", showWarnings = FALSE)
