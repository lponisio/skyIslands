library(bipartite, quietly = TRUE)
library(lme4, quietly = TRUE)
library(lmerTest, quietly = TRUE)
library(ggplot2, quietly = TRUE)
source('src/misc.R')
source('src/calcPca.R')
source('src/calcSpec.R')


load('../../data/nets.Rdata')
load('../../data/spec.Rdata')

save.path <- 'saved'
