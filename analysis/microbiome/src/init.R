library(ggplot2)
library(brms)
## library(tidybayes)
## library(tidyverse)
library(dplyr)


load('../../data/spec_RBCL_16s.Rdata')
site.sum <- read.csv("../../data/sitestats.csv")
veg <- read.csv("../../data/veg.csv")
load("../../data/trees.Rdata")

parasites <- c("AspergillusSpp",
               "AscosphaeraSpp",
               "ApicystisSpp",
               "CrithidiaExpoeki",
               "CrithidiaBombi",
               "NosemaBombi",
               "NosemaCeranae")

spec.net <- merge(spec.net, site.sum)

spec.net <- merge(spec.net, veg, all.x=TRUE)


dir.create(path="saved", showWarnings = FALSE)
dir.create(path="saved/tables", showWarnings = FALSE)
