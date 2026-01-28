library(bipartite, quietly = TRUE)
library(lme4, quietly = TRUE)
library(lmerTest, quietly = TRUE)
source('analysis/role/src/misc.R')
source('analysis/role/src/calcPca.R')
source('analysis/role/src/calcSpec.R')

load('data/spec_net.Rdata')
load('data/networks/Year_PlantPollinator_Bees.Rdata')
save.path <- 'saved'

save.dir.comm <- "saved/communities"
save.dir.nulls <- "saved/nulls"

type <- "plant"
nnull <- 9
species.type="PlantGenusSpecies" #changed from specimen (GenusSpecies)
species.type.int="GenusSpecies" #changed from plant (PlantGenusSpecies)


