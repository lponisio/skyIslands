
betalinkPP <- function (n1, n2, bf = B01, lower.level, higher.level) {
    ## From Poisot et
    ## al. 2012 https://doi.org/10.1111/ele.12002
    ## function modified from betalink packages that including
    ## calculating turnover of lower.level and pollinators seperatly
    v1 <- igraph::V(n1)$name
    v2 <- igraph::V(n2)$name
    vs <- v1[v1 %in% v2]
    beta_S <- bf(betapart(v1, v2))
    beta_S_lower.level <- bf(betapart(v1[v1 %in% lower.level], v2[v2 %in%
                                                                  lower.level]))
    beta_S_higher.level <- bf(betapart(v1[v1 %in% higher.level], v2[v2 %in%
                                                                    higher.level]))
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
                S_lower.level=beta_S_lower.level,
                S_higher.level=beta_S_higher.level,
                OS = beta_OS,
                WN = beta_WN,
                ST = beta_ST))
}


networkBetadiversity <- function (N, complete = FALSE,
                                  lower.level, higher.level,
                                  geo.dist,
                                  nets.by.SR=FALSE, ...) {
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
                            lower.level=lower.level, higher.level=higher.level, ...)
            b$i <- names(N)[i]
            b$j <- names(N)[j]
            beta <- rbind(beta, rbind(unlist(b)))
        }
    }
    if (NROW(beta) == 1) {
        beta <- data.frame(t(beta[, c("i", "j", "S", "OS", "WN",
                                      "ST",  "S_lower.level", "S_higher.level")]))
    }
    else {
        beta <- data.frame(beta[, c("i", "j", "S", "OS", "WN",
                                    "ST", "S_lower.level", "S_higher.level")])
    }
    beta$OS <- as.numeric(as.vector(beta$OS))
    beta$S <- as.numeric(as.vector(beta$S))
    beta$S_higher.level <- as.numeric(as.vector(beta$S_higher.level))
    beta$S_lower.level <- as.numeric(as.vector(beta$S_lower.level))
    beta$WN <- as.numeric(as.vector(beta$WN))
    beta$ST <- as.numeric(as.vector(beta$ST))
    beta$PropST <- beta$ST/beta$WN

    ## site and years names
    beta$Site1 <- sapply(strsplit(as.character(beta$i), "\\."),
                         function(x) x[[1]])
    beta$Site2 <- sapply(strsplit(as.character(beta$j), "\\."),
                         function(x) x[[1]])
    # beta$Year1 <- sapply(strsplit(as.character(beta$i), "\\."),
    #                      function(x) x[[2]])
    # beta$Year2 <- sapply(strsplit(as.character(beta$j), "\\."),
    #                      function(x) x[[2]])

    if(nets.by.SR){
        beta$SR1 <- sapply(strsplit(as.character(beta$i), "\\."),
                           function(x) x[[3]])
        beta$SR2 <- sapply(strsplit(as.character(beta$j), "\\."),
                           function(x) x[[3]])
        beta$Date1 <- spec$Date[match(paste0(beta$SR1, beta$Year1),
                                      paste0(spec$SampleRound, spec$Year))]
        beta$Date2 <- spec$Date[match(paste0(beta$SR2, beta$Year2),
                                      paste0(spec$SampleRound,
                                             spec$Year))]
        beta$SRDist <- apply(beta, 1, function(x){
            as.numeric(as.Date(x["Date2"]) -  as.Date(x["Date1"]))
        })
    }


    beta$GeoDist <- apply(beta, 1, function(x){
        geo.dist[x["Site1"],  x["Site2"]]
    })

    beta$YrDist <- apply(beta, 1, function(x){
        as.numeric(x["Year2"]) -  as.numeric(x["Year1"])
    })


    return(beta)
}

