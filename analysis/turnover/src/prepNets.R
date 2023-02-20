
## the purpose of this function is to break up data with many
## sites/years and prepare it for network analysis.

dropNet <- function(z){
  z[!sapply(z, FUN=function(q){
    any(dim(q) < 3)
  })]
}


###  adj matrices by site, yr, SR
makeNets <- function(spec.dat, net.type,
                     species=c("Plant", "Pollinator"), ...){
    spec.dat$YearSR <- paste(spec.dat$Year, spec.dat$SampleRound, sep=".")
    nets <- breakNet(spec.dat, 'Site', 'YearSR', ...)
    nets <- lapply(nets, bipartite::empty)

    ## graphs
    nets.graph <- lapply(nets, graph_from_incidence_matrix,
                         weighted =   TRUE, directed = FALSE)
    nets.graph <-  lapply(nets.graph, function(x){
        vertex_attr(x)$type[vertex_attr(x)$type] <- species[2]
        vertex_attr(x)$type[vertex_attr(x)$type
                            != species[2]] <- species[1]
        return(x)
    })

    ## unweighted
    nets.graph.uw <- lapply(nets, graph_from_incidence_matrix,
                            directed = FALSE)
    nets.graph.uw <-  lapply(nets.graph.uw, function(x){
        vertex_attr(x)$type[vertex_attr(x)$type] <- species[2]
        vertex_attr(x)$type[vertex_attr(x)$type
                            != species[2]] <- species[1]
        return(x)
    })

    years <- sapply(strsplit(names(nets), "[.]"), function(x) x[[2]])
    sites <- sapply(strsplit(names(nets), "[.]"), function(x) x[[1]])

    save(nets.graph,nets.graph.uw, nets, years, sites,
         file=sprintf("../data/nets%s%s.Rdata", net.type,
                      paste(species, collapse="")
                      ))

    ## species stats
    sp.lev <- calcSpec(nets)
    save(sp.lev, file=sprintf('../data/splev%s%s.Rdata', net.type,
                      paste(species, collapse="")
                      ))
}



breakNet <- function(spec.dat, site, year,
                     higher.level="GenusSpecies",
                     lower.level="PlantGenusSpecies",
                     mean.by.year=FALSE ){
    ## breaks network by site and year, and if mean.by.year=TRUE,
    ## takes the mean across sampling rounds
    ## puts data together in a list and removes empty matrices
    if(lower.level == "PlantGenusSpecies"){
    agg.spec <- aggregate(list(abund=spec.dat[, higher.level]),
                          list(HigherLevel=spec.dat[, higher.level],
                               Site=spec.dat[,site],
                               Year=spec.dat[,year],
                               LowerLevel=
                                   spec.dat[, lower.level]),
                          length)
    } else if(higher.level == "Parasite"){
          agg.spec <- aggregate(list(abund=spec.dat[, "count"]),
                          list(HigherLevel=spec.dat[, higher.level],
                               Site=spec.dat[,site],
                               Year=spec.dat[,year],
                               LowerLevel=
                                   spec.dat[, lower.level]),
                          mean, na.rm=TRUE)
    }

    if(mean.by.year){
        yrs <- sapply(strsplit(agg.spec$Year, "[.]"), function(x){x[[1]]})
        agg.spec <- aggregate(list(abund=agg.spec$abund),
                              list(HigherLevel=agg.spec$HigherLevel,
                                   Site=agg.spec$Site,
                                   Year=yrs,
                                   LowerLevel=
                                       agg.spec$LowerLevel),
                              mean)
    }
    sites <- split(agg.spec, agg.spec$Site)
    networks <- lapply(sites, function(x){
        split(x, f=x[,"Year"])
    })
    ## formats data matrices appropriate for network analysis
    comms <- lapply(unlist(networks, recursive=FALSE), function(y){
        samp2site.spp(site=y[,"LowerLevel"],
                      spp=y[,"HigherLevel"],
                      abund=y[,"abund"])
    })
    return(comms)
}


getSpecies <- function(networks, FUN){
  species.site <- lapply(networks, FUN)
  site.plant <- rep(names(species.site), lapply(species.site, length))
  species <- data.frame(species=do.call(c, species.site),
                        siteStatus=site.plant,
                        site= sapply(strsplit(site.plant, "_"),
                                     function(x) x[1]),
                        status= sapply(strsplit(site.plant, "_"),
                                       function(x) x[2]))
  return(species)
}



## calculates the degree of species in a network
getDegree <- function(x, MARGIN){
  apply(x, MARGIN, function(y){
    length(y[y != 0])/length(y)
  })
}

## calculate various stats
calcStats <- function(x){
  means=mean(x)
  medians=median(x)
  mins <- min(x)
  maxs <- max(x)
  sds <- sd(x)
  return(c(mean=means,
           median=medians,
           min=mins,
           max=maxs,
           sd=sds))
}



## number of species that interact
getCon <- function(x, INDEX){
  apply(x, INDEX, function(y) sum(y > 0))
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
  all.pp$speciesType <- c(rep("pollinator", n.pp[1]),
                          rep("plant", n.pp[2]))
  all.pp$Site <- strsplit(names.net, seps)[[1]][1]
  all.pp$assem <- strsplit(names.net, seps)[[1]][2]
  return(all.pp)
}

