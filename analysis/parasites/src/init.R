library(ggplot2)
library(brms)
## library(tidybayes)
library(tidyverse)


load('../../data/spec_RBCL_16s.Rdata')
site.sum <- read.csv("../../data/sitestats.csv")
veg <- read.csv("../../data/veg.csv")


parasites <- c("AspergillusSpp",
               "AscosphaeraSpp",
               "ApicystisSpp",
               "CrithidiaExpoeki",
               "CrithidiaBombi",
               "NosemaBombi",
               "NosemaCeranae")

## subset to the bees we screened and the screenings worked
## spec <- spec[spec$Apidae == 1 &
##              !is.na(spec$Apidae),]

spec <- merge(spec, site.sum)

spec <- merge(spec, veg, all.x=TRUE)
