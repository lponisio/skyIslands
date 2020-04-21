library(igraph, quietly = TRUE)
library(bipartite, quietly = TRUE)
library(SYNCSA, quietly = TRUE)

calcMetric <- function(dat.web, ...) {
    mets <-  networklevel(dat.web, ...)
    ## the functional redundancy function takes a matrix of sites and
    ## species, and a trait matrix whwere the rownames of the traits
    ## match the column names of the site by species matric. In our
    ## case, the plants and pollinators are the "species"
    ## respectively, and their traits are their interaction partners.

    ## create names for later use in site x species matrix
    rownames(dat.web) <- 1:nrow(dat.web)
    colnames(dat.web) <- 1:ncol(dat.web)

    ## site by species matrix where there is only one "site" since the
    ## networks are site specific, and the columns are the
    ## species.

    ## abundance weighted site by species matrices
    plants <- matrix(rowSums(dat.web),  nrow=1)
    pols <- matrix(colSums(dat.web),  nrow=1)

    colnames(plants) <- rownames(dat.web)
    colnames(pols) <- colnames(dat.web)
    rownames(plants) <- "Plants"
    rownames(pols) <- "Pols"

    ## pull out Functional redundancy score based on: de Bello, F.;
    ## Leps, J.; Lavorel, S. & Moretti, M. (2007). Importance of
    ## species abundance for assessment of trait composition: an
    ## example based on pollinator communities. Community Ecology, 8,
    ## 163:170. and functional complementarity score based on Rao,
    ## C.R. (1982). Diversity and dissimilarity coefficients: a
    ## unified approach. Theoretical Population Biology, 21, 24:43.

    redund.plant <- unlist(rao.diversity(plants,
                                         traits=
                                             dat.web)[c("FunRao",
                                                        "FunRedundancy")])
    redund.pol <- unlist(rao.diversity(pols,
                                       traits=
                                           t(dat.web))[c("FunRao",
                                                         "FunRedundancy")])

    return(c(mets,
             redund.plant,
             redund.pol))

}


## function to simulate 1 null, and calculate statistics on it
calcNullStat <- function(dat.web,
                         null.fun,...) {
    sim.web <- null.fun(dat.web)
    out.met <- calcMetric(sim.web,...)
    return(out.met)
}

##  function that computes summary statistics on simulated null matrices
##  (nulls simulated from web N times)
calcNetworkMetrics <- function(dat.web, N,
                               index= c("niche overlap",
                                        "functional complementarity",
                                        "weighted NODF",
                                        "ISA",
                                        "SA",
                                        "vulnerability",
                                        "number of species",
                                        "cluster coefficient",
                                        "H2")) {
    ## calculate pvalues
    pvals <- function(stats, nnull){
        rowSums(stats >= stats[, rep(1, ncol(stats))])/(nnull + 1)
    }
    ## calculate zvalues two different ways
    zvals <-function(stats){
        z.sd <- (stats[,1] -
                 apply(stats, 1, mean, na.rm = TRUE))/
            apply(stats, 1, sd, na.rm = TRUE)
        z.sd[is.infinite(z.sd)] <- NA
        return(z.sd)
    }
    ## check that matrix is proper format (no empty row/col and no NAs)
    ## drop empty rows and columns
    dat.web <- as.matrix(bipartite::empty(dat.web))
    ## check to make sure emptied matrix is large enough
    ## to calculate statistics on
    if(all(dim(dat.web) >= 2)) {
        ## calculate null metrics
        null.stat <- replicate(N,
                               calcNullStat(dat.web,
                                            null.fun= vaznull.fast,
                                            index=index,
                                            dist="chao"),
                               simplify=TRUE)
        ## calculate metrics from data
        true.stat <- calcMetric(dat.web,
                                index=index)
        out.mets <- cbind(true.stat, null.stat)
        ## compute z scores
        zvalues <- zvals(out.mets)
        names(zvalues) <- paste("z", names(true.stat), sep="")
        ## compute p-values
        pvalues <- pvals(out.mets, N)
        names(pvalues) <- paste("p", names(true.stat), sep="")
        out <- c(true.stat, zvalues, pvalues)
        return(out)


    }
    return(rep(NA, (length(index) + 6)*3))
}

prepDat <- function(cor.stats, spec.dat,
                    cols.to.keep= c("Site", "Lat", "LatPoly1",
                                    "LatPoly2", "Elev", "Area"), net.type){
    dats <- do.call(rbind, cor.stats)
    out <- data.frame(dats)
    out$Site <- sapply(strsplit(names(cor.stats), "\\."),
                       function(x) x[1])
    out$Year <-  as.factor(sapply(strsplit(names(cor.stats), "\\."),
                                  function(x) x[2]))

    if(net.type == "YrSR"){
        out$SampleRound <-  as.factor(sapply(strsplit(names(cor.stats), "\\."),
                                             function(x) x[3]))
        temporal.dat <- unique(spec.dat[, c("Site","SampleRound", "Year",
                                            "DoyPoly1",
                                            "DoyPoly2", "Date", "Doy")])
        out <- merge(out, temporal.dat)
    }
    site.dats <- unique(spec[, c(cols.to.keep)])
    site.dats$Site  <- as.character(site.dats$Site)
    site.dats <- apply(site.dats[, cols.to.keep[-1]], 2, function(x)
        tapply(x, site.dats$Site, mean, na.rm=TRUE))
    site.dats <- as.data.frame(site.dats)
    site.dats$Site <- rownames(site.dats)
    out <- merge(out, site.dats, by="Site", all.x=TRUE)
    out$Site <- factor(out$Site, levels=c("CH", "PL",
                                          "HM", "MM", "SC", "SM", "JC"))
    rownames(out) <- NULL
    return(out)
}
