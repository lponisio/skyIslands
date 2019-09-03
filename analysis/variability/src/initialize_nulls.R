library(bipartite, quietly = TRUE)
source('src/misc.R')
source('src/vaznull2.R')
source('src/commPrep.R')

load('../../data/spec.Rdata')

save.dir.comm <- "saved/communities"
save.dir.nulls <- "saved/nulls"


sites <- unique(spec$Site)
years <- unique(spec$Year)
spec$Int <- paste(spec$GenusSpecies,
                  spec$PlantGenusSpecies)

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 0){
    sp.type <- args[1]
    net.type <- args[2]
    nnull <- args[3]
} else{
    sp.type <- "pol"
    net.type <- "Site"
    nnull <- 99
}

if(sp.type=="pol"){
    species.type="GenusSpecies"
    species.type.int="PlantGenusSpecies"
}

if(sp.type=="plants"){
    species.type="PlantGenusSpecies"
    species.type.int="GenusSpecies"
}

if(net.type == "Year"){
    to.lapply <- years
}

if(net.type == "Site"){
    to.lapply <- sites
}
