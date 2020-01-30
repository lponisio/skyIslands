## calculates the species roles from a network and returns a dataframe
## with site status and ypr
## takes networks and specimen data

calcSpec <- function(nets, dist.metric="gower"){
    ## applies specieslevel from bipartite to networks
    nets <- nets[sapply(nets, function(x){
        all(dim(x) > 2)
    })]
    species.lev <- lapply(nets, function(x){
        sl <- specieslevel(x)
        sl$'higher level'$tot.int <- colSums(x)[colSums(x) != 0]
        sl$'lower level'$tot.int <- rowSums(x)[rowSums(x) != 0]
        lower.level.niche.overlap <- as.matrix(vegdist(x,
                                                       method=dist.metric))
        diag(lower.level.niche.overlap) <- NA
        sl$'lower level'$niche.overlap <- apply(lower.level.niche.overlap, 1, mean,
                                                na.rm=TRUE)

        higher.level.niche.overlap <- as.matrix(vegdist(t(x),
                                                        method=dist.metric))
        diag(higher.level.niche.overlap) <- NA
        sl$'higher level'$niche.overlap <- apply(higher.level.niche.overlap, 1, mean,
                                                 na.rm=TRUE)
        sl$'higher level'$rare.degree <- apply(x, 2, chao1)
        sl$'lower level'$rare.degree <- apply(x, 1, chao1)
        return(sl)
    })

    ## extract the values and make a dataframe
    specs  <-  mapply(function(a, b)
        getSpec(species.lev = a,
                names.net = b,
                seps="[.]"),
        a = species.lev,
        b = names(nets),
        SIMPLIFY = FALSE)

    specs <- do.call(rbind, specs)
    rownames(specs) <- NULL
    return(specs)
}


## extreact specialization scores from specieslevel function and
## return data frame
getSpec <- function(species.lev, names.net, seps="_"){
    n.pp <- sapply(species.lev, nrow)
    pp <- c(unlist(sapply(species.lev, rownames)))
    names(pp) <- NULL
    all.pp <- do.call(rbind, species.lev)
    rownames(all.pp) <- NULL
    try(all.pp$GenusSpecies <- pp)
    all.pp$speciesType <- c(rep("higher.level", n.pp[1]),
                            rep("lower.level", n.pp[2]))
    all.pp$Site <- strsplit(names.net, seps)[[1]][1]
    all.pp$Year <- strsplit(names.net, seps)[[1]][2]
    all.pp$SpSiteYear <- paste0(all.pp$GenusSpecies, all.pp$Site, all.pp$Year)

    SR <- try(strsplit(names.net, seps)[[1]][3], silent=TRUE)
    if(!is.na(SR)){
        all.pp$SR <- SR
        all.pp$SpSiteYear <- paste0(all.pp$GenusSpecies, all.pp$Site,
                                    all.pp$Year, all.pp$SR)
    }
    return(all.pp)
}

