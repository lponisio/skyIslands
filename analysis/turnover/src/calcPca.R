calcNetworkPca <- function(species.roles,
                           loadings,
                           metrics,
                           nets.by.SR=FALSE){
    ## this function calculates a pca of network roles, then takes the
    ## mean, and var (methods passed in as arguments) for each
    ## species. It returns the mean and variable for each species, as
    ## well as the pca loadings
    species.roles <-  species.roles[!apply( species.roles[,metrics], 1,
                                           function(x) any(is.na(x))),]
    mets.only <-  species.roles[, metrics]
    ## runs the pca
    all.pca <- prcomp(mets.only, scale. = TRUE, center = TRUE)
    ## make a nice dataframe
    if(nets.by.SR){
        all.spp <-  species.roles[, c("Site", "Year",
                                      "GenusSpecies", "SR",
                                      "SpSiteYear")]
    } else{
        all.spp <-  species.roles[, c("Site", "Year",
                                      "GenusSpecies",
                                      "SpSiteYear")]
    }

    all.spp <- cbind(all.spp,
                     pca=all.pca$x[,loadings])

    return(list(pcas=all.spp,
                pca.loadings = all.pca))
}

calcDiffPcas <-  function(x, geo.dist, nets.by.SR=FALSE){
    ## function to calculate all PCA scores differences within a
    ## species
    ## all combinations of pcas
    if(dim(x)[1] > 1){
        pcas <- expand.grid(pca1=x$pca, pca2=x$pca)
        ## and their site year combos
        if(nets.by.SR){
            pcas <- cbind(pcas,
                          expand.grid(SiteYear1=paste(x$Site,
                                                      x$Year, x$SR, sep=":"),
                                      SiteYear2=paste(x$Site,
                                                      x$Year, x$SR,
                                                      sep=":")))
        } else{
            pcas <- cbind(pcas,
                          expand.grid(SiteYear1=paste(x$Site,
                                                      x$Year,  sep=":"),
                                      SiteYear2=paste(x$Site,
                                                      x$Year,
                                                      sep=":")))
        }
        ## from the diagonal (same pca comarison)
        pcas <- pcas[pcas$SiteYear1 != pcas$SiteYear2,]
        ## take the difference
        pcas$diffPca <- abs(pcas$pca1 - pcas$pca2)

        ## clean up output
        pcas$Site1 <- sapply(strsplit(as.character(pcas$SiteYear1), ":"),
                             function(x) x[[1]])
        pcas$Year1 <- sapply(strsplit(as.character(pcas$SiteYear1), ":"),
                             function(x) x[[2]])

        pcas$Site2 <- sapply(strsplit(as.character(pcas$SiteYear2), ":"),
                             function(x) x[[1]])
        pcas$Year2 <- sapply(strsplit(as.character(pcas$SiteYear2), ":"),
                             function(x) x[[2]])

        if(nets.by.SR){
            pcas$SR1 <- sapply(strsplit(as.character(pcas$SiteYear1), ":"),
                               function(x) x[[3]])
            pcas$SR2 <- sapply(strsplit(as.character(pcas$SiteYear2), ":"),
                               function(x) x[[3]])
        }

        pcas$SiteYear1 <- NULL
        pcas$SiteYear2 <- NULL

        pcas$GeoDist <- apply(pcas, 1, function(x){
            geo.dist[x["Site1"],  x["Site2"]]
        })
        ## add back species label
        print(unique(x$GenusSpecies))
        ## if(unique(x$GenusSpecies) == "Calliopsis sp. a") browser()
        pcas$GenusSpecies <- unique(x$GenusSpecies)
    } else {
        pcas <- NULL
    }
    return(pcas)
}
