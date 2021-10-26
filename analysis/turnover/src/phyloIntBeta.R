

calcCommDis <- function(spec.dat, type, lc, abund.w = TRUE,
                        between, within,
                        geo.dist = NA ){
    ## prep community of sites and interactions phylo-beta of
    ## interactions takes raw specimen data (spec.dat), the community
    ## type, the linkcommunity object from linkcomm, and whether to
    ## weight by abundance

    ## within  is the column to group species by. If == "year" then
    ## within a site across years, if == "site" when within a
    ## year between sites

    ## dist.site is the geographic distances between sites if
    ## computing between sites within a year, NA else
browser()
    prep.comm <- aggregate(spec.dat[, type],
                           list(site=spec.dat[,within],
                                year=spec.dat[,between],
                                sp=spec.dat[, type]),
                           length)
    ## split  by site
    bysite <- split(prep.comm, prep.comm$site)
    ## create a site by interaction community
    comm <- lapply(bysite, function(y) {
        samp2site.spp(y$year, y$sp, y$x)
    })
    ## drop sites that were only samples one year
    comm <- comm[lapply(comm, nrow) >= 2]

    ## calculate the number of years between sampling times
    if(between == "Year"){
        years <- lapply(comm, rownames)
        dist.time <- as.matrix(dist(unique(spec.dat$Year), diag=TRUE))
        rownames(dist.time) <- colnames(dist.time) <- unique(spec.dat$Year)
        dists.yrs <- lapply(years, function(x){dist.time[x,x]})
        c.dist <- lapply(dists.yrs, function(x)
        {as.matrix(x)[lower.tri(x)]})
     }else{
        sites <- lapply(comm, rownames)[[1]]
        geo.dist <- geo.dist[sites, sites]
        c.dist <- geo.dist[lower.tri(geo.dist)]
     }
    ## rename the species as node numbers in preparation for using
    ## comdist using the dendrogram as distances

    node.num <- cbind(lc$hclust$labels, paste(lc$edgelist[,1],
                                              lc$edgelist[,2]))
    node.order <- lapply(comm, function(x){
        node.num[,1][match(colnames(x), node.num[,2])]
    })
    assignCols <- function(a, b){
        colnames(a) <- b
        return(a)
    }
    comm.int <- mapply(function(a, b)
        assignCols(a,b),
        a=comm,
        b=node.order)

    ## calculate "phylogenetic" distance between interactions
    dist.tree <- cophenetic(lc$hclust)

    ## calculate the dissimilarity between years based on the dengrogram
    comm.dis <- lapply(comm.int, comdist,
                       dis=dist.tree,
                       abundance.weighted= abund.w)

    ## return a dataframe
    c.comm <- lapply(comm.dis, function(x) {
        as.matrix(x)[lower.tri(x)]
    })
    site <- names(comm.dis)
    sites <- rep(site, sapply(c.comm, length))
    comp <- lapply(comm.dis, function(x){
        x <- as.matrix(x)
        comb <- combn(rownames(x), 2)
        out <- apply(comb, 2, function(x) paste(unlist(x), collapse =" "))
        return(out)
    })

    out <- as.data.frame(cbind(PhyloInt=unlist(c.comm),
                               Dist= unlist(c.dist)))
    out$Site <- as.factor(sites)
    out$Comp <- unlist(comp)
    rownames(out) <- NULL
    return(list(phylo.int = out,
                dist.tree=dist.tree,
                comm=comm))
}


calcDendDis <- function(spec.dat, type){
    ## created a dendrograph of shared interactions between nodes
    prep.comm <- aggregate(spec.dat[, type[1]],
                           list(sp=spec.dat[, type[1]],
                                sp2= spec.dat[, type[2]]), length)
    comm <- samp2site.spp(prep.comm$sp, prep.comm$sp2, prep.comm$x)
    comm.dis <- vegdist(comm, "chao", diag= TRUE)
    dengram <- hclust(comm.dis, method="average")

    fden <- function(){
        plot(dengram, cex.lab=0.5, hang=-1)
    }
    path.dir <- 'figures'
    pdf.f(fden, file= file.path(path.dir, sprintf("%s.pdf",
                                                  paste("dedrogram",
                                                        type[1],
                                                        sep="_"))),
          width=25, height=10)
    return(dengram)
}
