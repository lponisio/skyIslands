args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 0){
    sp.type <- args[1]
    net.type <- args[2]
    occ <- args[3]
} else{
    sp.type <- "pol"
    net.type <- "Site"
    occ <- "abund"
}

library(vegan)
source('src/misc.R')
source('src/beta.R')

if(sp.type == "pol"){
    speciesType <- "pollinator"
} else{
    speciesType <- "plants"
}


if(occ == "abund"){
    binary <- FALSE
    dis.method <- "gower"
    load(file=file.path('saved/communities',
                        sprintf('%s-%s-abund.Rdata',
                                sp.type,
                                net.type)))
    load(file=file.path('saved/nulls',
                        sprintf('%s-%s-alpha.Rdata',
                                sp.type,
                                net.type)))
}

if(occ == "occ"){
    occ <- "occ"
    binary <- TRUE
    dis.method <- "jaccard"
    load(file=file.path('saved/communities',
                        sprintf('%s-%s-abund.Rdata',
                                sp.type,
                                net.type)))
    load(file=file.path('saved/nulls',
                        sprintf('%s-%s-occ.Rdata',
                                sp.type,
                                net.type)))
}

if(sp.type=="pol"){
    ylabel <- "Pollinator species turnover"
}
if(sp.type=="ints"){
    ylabel <- "Interaction turnover"
}
if(sp.type=="plants"){
    ylabel <- "Plant species turnover"
}
