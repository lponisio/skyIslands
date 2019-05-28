betalink_test  <- function (n1, n2, bf) {
    browser()
    v1 <- igraph::V(n1)$name
    v2 <- igraph::V(n2)$name
    vs <- v1[v1 %in% v2]
    beta_S <- bf(betapart(v1, v2))
    e1 <- plyr::aaply(igraph::get.edgelist(n1), 1, function(x) stringr::str_c(x,
        collapse = "--", paste = "_"))
    e2 <- plyr::aaply(igraph::get.edgelist(n2), 1, function(x) stringr::str_c(x,
        collapse = "--", paste = "_"))
    beta_WN <- bf(betapart(e1, e2))
    if (length(vs) >= 2) {
        sn1 <- igraph::induced.subgraph(n1, which(igraph::V(n1)$name %in%
            vs))
        sn2 <- igraph::induced.subgraph(n2, which(igraph::V(n2)$name %in%
            vs))
        se1 <- plyr::aaply(igraph::get.edgelist(sn1), 1, function(x) stringr::str_c(x,
            collapse = "--", paste = "_"))
        se2 <- plyr::aaply(igraph::get.edgelist(sn2), 1, function(x) stringr::str_c(x,
            collapse = "--", paste = "_"))
        beta_OS <- bf(betapart(se1, se2))
        beta_ST <- beta_WN - beta_OS
    }
    else {
        beta_OS <- NaN
        beta_ST <- NaN
    }
    return(list(S = beta_S, OS = beta_OS, WN = beta_WN, ST = beta_ST))
}


network_betadiversity_test <- function (N, complete = FALSE, ...) {
    N <- name_networks(N)
    beta <- NULL
    stop_at <- ifelse(complete, length(N), length(N) - 1)
    for (i in c(1:stop_at)) {
        start_at <- ifelse(complete, 1, i + 1)
        inner_stop <- length(N)
        inner_steps <- c(start_at:inner_stop)
        inner_steps <- inner_steps[which(inner_steps != i)]
        for (j in inner_steps) {
            b <- betalink_test(N[[i]], N[[j]], ...)
            b$i <- names(N)[i]
            b$j <- names(N)[j]
            beta <- rbind(beta, rbind(unlist(b)))
        }
    }
    if (NROW(beta) == 1) {
        beta <- data.frame(t(beta[, c("i", "j", "S", "OS", "WN",
            "ST")]))
    }
    else {
        beta <- data.frame(beta[, c("i", "j", "S", "OS", "WN",
            "ST")])
    }
    beta$OS <- as.numeric(as.vector(beta$OS))
    beta$S <- as.numeric(as.vector(beta$S))
    beta$WN <- as.numeric(as.vector(beta$WN))
    beta$ST <- as.numeric(as.vector(beta$ST))
    return(beta)
}

