args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 0){
    sp.type <- args[1]
    net.type <- args[2]
    nnull <- args[3]
} else{
    sp.type <- "pol"
    net.type <- "Year"
    nnull <- 999
}


library(bipartite, quietly = TRUE)
source('src/misc.R')
source('src/probNull.R')
source('src/commPrep.R')

load('../../data/spec_net.Rdata')
spec <- spec.net

save.dir.comm <- "saved/communities"
save.dir.nulls <- "saved/nulls"

spec$SiteYear <- paste(spec$Site, spec$Year, sep=":")
sites <- unique(spec$Site)
years <- unique(spec$Year)
spec$Int <- paste(spec$GenusSpecies,
                  spec$PlantGenusSpecies)
site.year <- unique(spec$SiteYear)

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

if(net.type == "SiteYear"){
    to.lapply <- site.year
}

