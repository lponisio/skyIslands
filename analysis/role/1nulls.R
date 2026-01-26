## setwd('~/Dropbox/skyislands')
rm(list=ls())
setwd('role')
source('src/initialize_nulls.R')

## ************************************************************
## create community matrices for species species by year or site
## ************************************************************

## net.type is either, Site, Year or SiteYear for looking at turnover
## between years withina  a site, between sites within a year, and ???

## when lapplying on years, generates networks for each species
## (pollinator if species.type=="GenusSpecies) to examine turnover
## between sites within a year

## when lapplying on sites, generates networks for each species
## (pollinator if species.type=="GenusSpecies) to examine turnover
## between years within a site

comms <- lapply(to.lapply, calcSiteBeta,
                species.type=species.type,
                spec=spec,
                species.type.int=species.type.int,
                type.net= net.type)

comm <- makePretty(comms, spec, net.type)

save(comm, file=file.path(save.dir.comm,
                          sprintf('%s-%s-abund.Rdata',
                                  sp.type,
                                  net.type)))

## ************************************************************
## alpha div nulls
## ************************************************************
nulls <- rapply(comm$comm, calcProbNull, N=nnull, how="replace")

save(nulls, file=file.path(save.dir.nulls,
                           sprintf('%s-%s-alpha.Rdata',
                                   sp.type,
                                   net.type)))

## ************************************************************
## occurrence nulls
## ************************************************************
occ.null <- function(web){
    simulate(vegan::nullmodel(web, method="quasiswap"),1)[,,1]
}

rep.occ.null <- function(web, N){
    replicate(N, occ.null(web), simplify = FALSE)
}

nulls <- rapply(comm$comm, rep.occ.null, N=nnull, how="replace")

save(nulls, file=file.path(save.dir.nulls,
                           sprintf('%s-%s-occ.Rdata',
                                   sp.type,
                                   net.type)))

## ************************************************************
## test
## ************************************************************
## set interactions to be all the same across sampling rounds for each
## species
if(net.type == "Year"){
    test.comm <- lapply(comm$comm$'2012', function(x){
        x[1:nrow(x), 1:ncol(x)] <- 0
        x[,c(1,2,3)] <- 2
        x
    })
    comm$comm$'2012' <- test.comm
    nulls <- rapply(comm$comm, calcProbNull, N=nnull, how="replace")

    save(nulls, file=file.path(save.dir.nulls,
                               sprintf('%s-%s-alpha-test.Rdata',
                                       sp.type,
                                       net.type)))
    save(comm, file=file.path(save.dir.comm,
                              sprintf('%s-%s-abund-test.Rdata',
                                      sp.type,
                                      net.type)))
}
