library(ggplot2)
library(brms)
## library(tidybayes)
library(tidyverse)


load('../../data/spec_RBCL_16s.Rdata')
site.sum <- read.csv("../../data/sitestats.csv")
veg <- read.csv("../../data/veg.csv")
load("../../data/trees.RData")

parasites <- c("AspergillusSpp",
               "AscosphaeraSpp",
               "ApicystisSpp",
               "CrithidiaExpoeki",
               "CrithidiaBombi",
               "NosemaBombi",
               "NosemaCeranae")

spec <- merge(spec, site.sum)

spec <- merge(spec, veg, all.x=TRUE)


dir.create(path="saved", showWarnings = FALSE)
dir.create(path="saved/tables", showWarnings = FALSE)
