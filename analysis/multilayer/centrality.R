## setwd('~/Dropbox/skyIslands')
rm(list=ls())
setwd('analysis/multilayer')

library(multinet)
library(ggplot2)
library(bipartite)
library(igraph)

load('../../data/spec.Rdata')
load('../../data/netsYr.Rdata')
load('../../data/splevYr.Rdata')


ml.nets.files <-  list.files("saved", "mpx")

ml.nets <- lapply( file.path("saved", ml.nets.files), read_ml)
names(ml.nets) <- sapply(strsplit(ml.nets.files, "[.]"), function(x) x[[1]])

## net.2012.adj <-  as.matrix(as_adjacency_matrix(net.2012.all))
## nets.2012.mets <- specieslevel(net.2012.adj)

getMlDegree <- function(net){
    deg.ml <- degree_ml(net)
    ml.degrees <- data.frame(MLdegree= deg.ml[order(-deg.ml)],
                             GenusSpecies=actors_ml(net)[order(-deg.ml)])
    return(ml.degrees)
}


ml.degrees <- lapply(ml.nets, getMlDegree)
all.ml.degrees <- do.call(rbind, ml.degrees)

all.ml.degrees$Year <- rep(names(ml.degrees),
                           times=sapply(ml.degrees, nrow))

all.ml.degrees$Year <- gsub("net", "", all.ml.degrees$Year)

sp.lev$mldegree <- all.ml.degrees$MLdegree[
                   match(paste(sp.lev$GenusSpecies, sp.lev$Year),
                   paste(all.ml.degrees$GenusSpecies,
                         all.ml.degrees$Year))]


yvars <- c("S", "WN", "S_Plants", "S_Pols", "PropST", "OS")

ylabs <- c("Species turnover", "Interaction Turnover",
           "Plant turnover",
           "Pollinator turnover",
           "Interaction turnover due to species turnover",
           "Interaction turnover due to rewiring")

panels <- vector("list", length(yvars))

for(i in 1:length(yvars)){
    this.beta <- beta.same.year
    colnames(this.beta)[colnames(this.beta) == yvars[i]] <- "y"
    panels[[i]] <- ggplot(this.beta,
       aes(x=GeoDist, y=y,
           color=Year1)) + geom_point() + geom_smooth() +
    labs(x="Geographic Distance", y=ylabs[i]) +
    lims(y=c(0,1))
}

do.call(grid.arrange, panels)


