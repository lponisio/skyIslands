
betalinkPP <- function (n1, n2, bf = B01, plants, pols) {
    ## From Poisot et
    ## al. 2012 https://doi.org/10.1111/ele.12002
    ## function modified from betalink packages that including
    ## calculating turnover of plants and pollinators seperatly
    v1 <- igraph::V(n1)$name
    v2 <- igraph::V(n2)$name
    vs <- v1[v1 %in% v2]
    beta_S <- bf(betapart(v1, v2))
    beta_S_plant <- bf(betapart(v1[v1 %in% plants], v2[v2 %in%
                                                       plants]))
    beta_S_pol <- bf(betapart(v1[v1 %in% pols], v2[v2 %in%
                                                   pols]))
    e1 <- plyr::aaply(igraph::get.edgelist(n1), 1,
                      function(x) stringr::str_c(x,
                                                 collapse = "--", paste = "_"))
    e2 <- plyr::aaply(igraph::get.edgelist(n2), 1,
                      function(x) stringr::str_c(x,
                                                 collapse = "--", paste = "_"))
    beta_WN <- bf(betapart(e1, e2))
    if (length(vs) >= 2) {
        sn1 <- igraph::induced.subgraph(n1, which(igraph::V(n1)$name %in%
                                                               vs))
        sn2 <- igraph::induced.subgraph(n2, which(igraph::V(n2)$name %in%
                                                               vs))
        se1 <- plyr::aaply(igraph::get.edgelist(sn1), 1,
                           function(x) stringr::str_c(x,
                                                      collapse = "--", paste = "_"))
        se2 <- plyr::aaply(igraph::get.edgelist(sn2), 1,
                           function(x) stringr::str_c(x,
                                                      collapse = "--", paste = "_"))
        beta_OS <- bf(betapart(se1, se2))
        beta_ST <- beta_WN - beta_OS
    }
    else {
        beta_OS <- NaN
        beta_ST <- NaN
        beta_RW <- NaN
    }
    return(list(S = beta_S,
                S_Plants=beta_S_plant,
                S_Pols=beta_S_pol,
                OS = beta_OS,
                WN = beta_WN,
                ST = beta_ST))
}


network_betadiversity <- function (N, complete = FALSE,
                                   plants, pols, geo.dist,...) {
    ## modified function from betalink package. From Poisot et
    ## al. 2012 https://doi.org/10.1111/ele.12002

    ## this function mostly clean up the output of betalink, and add
    ## some data columns fo interest

    ## Bs: Dissimilarity in the species composition of communities
    ## Bwn: Dissimilarity of interactions
    ## Bos: Dissimilarity of interactions established between species
    ## common to both realisations
    ## Bst: Dissimilarity of interactions due to species turnover
    ## Bst/wn: Dissimilarity of interactions due to species turnover
    ## Bos/wn: propotion dissimilarity of interactions due to rewiring
    N <- name_networks(N)
    beta <- NULL
    stop_at <- ifelse(complete, length(N), length(N) - 1)
    for (i in c(1:stop_at)) {
        start_at <- ifelse(complete, 1, i + 1)
        inner_stop <- length(N)
        inner_steps <- c(start_at:inner_stop)
        inner_steps <- inner_steps[which(inner_steps != i)]
        for (j in inner_steps) {
            b <- betalinkPP(N[[i]], N[[j]],
                            plants=plants, pols=pols, ...)
            b$i <- names(N)[i]
            b$j <- names(N)[j]
            beta <- rbind(beta, rbind(unlist(b)))
        }
    }
    if (NROW(beta) == 1) {
        beta <- data.frame(t(beta[, c("i", "j", "S", "OS", "WN",
                                      "ST",  "S_Plants", "S_Pols")]))
    }
    else {
        beta <- data.frame(beta[, c("i", "j", "S", "OS", "WN",
                                    "ST", "S_Plants", "S_Pols")])
    }
    beta$OS <- as.numeric(as.vector(beta$OS))
    beta$S <- as.numeric(as.vector(beta$S))
    beta$S_Pols <- as.numeric(as.vector(beta$S_Pol))
    beta$S_Plants <- as.numeric(as.vector(beta$S_Plant))
    beta$WN <- as.numeric(as.vector(beta$WN))
    beta$ST <- as.numeric(as.vector(beta$ST))
    beta$PropST <- beta$ST/beta$WN

    ## site and years names
    beta$Site1 <- sapply(strsplit(as.character(beta$i), "\\."),
                         function(x) x[[1]])
    beta$Year1 <- sapply(strsplit(as.character(beta$i), "\\."),
                         function(x) x[[2]])
    beta$SR1 <- sapply(strsplit(as.character(beta$i), "\\."),
                       function(x) x[[3]])

    beta$Site2 <- sapply(strsplit(as.character(beta$j), "\\."),
                         function(x) x[[1]])
    beta$Year2 <- sapply(strsplit(as.character(beta$j), "\\."),
                         function(x) x[[2]])
    beta$SR2 <- sapply(strsplit(as.character(beta$j), "\\."),
                       function(x) x[[3]])

    beta$GeoDist <- apply(beta, 1, function(x){
        geo.dist[x["Site1"],  x["Site2"]]
    })

    beta$YrDist <- apply(beta, 1, function(x){
        as.numeric(x["Year2"]) -  as.numeric(x["Year1"])
    })


    return(beta)
}

