## calculates the species roles from a network and returns a dataframe
## with site status and ypr
## takes networks and specimen data  

calcSpec <- function(nets, spec){
  ## applies specieslevel from bipartite to networks
  species.lev <- lapply(nets, function(x){
    sl <- specieslevel(x)
    sl$'higher level'$tot.int <- colSums(x)
    sl$'lower level'$tot.int <- rowSums(x)
    sl$'higher level'$mean.k <- mean(sl$'higher level'$degree)
    sl$'lower level'$mean.k <- mean(sl$'lower level'$degree)
    sl$'higher level'$sd.k <- sd(sl$'higher level'$degree)
    sl$'lower level'$sd.k <- sd(sl$'lower level'$degree)
    sl$'higher level'$k <- (sl$'higher level'$degree -
                            sl$'higher level'$mean.k)/sl$'higher level'$sd.k
    sl$'lower level'$k <- (sl$'lower level'$degree -
                            sl$'lower level'$mean.k)/sl$'lower level'$sd.k
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
  all.pp$speciesType <- c(rep("pollinator", n.pp[1]),
                          rep("plant", n.pp[2]))
  all.pp$Site <- strsplit(names.net, seps)[[1]][1]
  all.pp$Year <- strsplit(names.net, seps)[[1]][2]
  return(all.pp)
}

