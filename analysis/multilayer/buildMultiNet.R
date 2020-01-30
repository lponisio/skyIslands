## setwd('~/Dropbox/skyIslands')
rm(list=ls())
setwd('analysis/multilayer')

library(multinet)
source("src/misc.R")
source("src/plotting.R")

net.type <- "Yr"
## species <- c("Plant", "Pollinator")
species <- c("Pollinator", "Parasite")

source('../turnover/src/initialize.R')

same.year <- split(nets.graph.uw, years)

spatial.order <- list('2012'=c("JC", "SC", "MM", "PL", "CH"),
                 '2017'= c("JC", "SC", "SM", "MM", "HM", "CH"),
                 '2018'= c("JC", "SC", "SM", "MM", "HM", "PL", "CH"))

yrs <- unique(years)

for(yr in yrs){
    print(yr)
    ml.net <- ml_empty()
    site.yr <- paste(spatial.order[[yr]], yr, sep=".")
    for(i in 1:length(spatial.order[[yr]])){
        add_igraph_layer_ml(ml.net, same.year[[yr]][[site.yr[i]]],
                            spatial.order[[yr]][i])
    }

    for(i in 1:(length(spatial.order[[yr]]) -1)){
        l1.to.l2 <- names(V(same.year[[yr]][[site.yr[i]]]))[
            names(V(same.year[[yr]][[site.yr[i]]])) %in%
            names(V(same.year[[yr]][[site.yr[i + 1]]]))]

        inter_layer_edges12 <- data.frame(
            actor1 = l1.to.l2,
            layer1 = rep(spatial.order[[yr]][i], length(l1.to.l2)),
            actor2 =  l1.to.l2,
            layer2 = rep(spatial.order[[yr]][i+1], length(l1.to.l2)))

        add_edges_ml(ml.net, inter_layer_edges12)

    }

    add_attributes_ml(ml.net, attributes="SpType", type="string",
                      target="actor")

    set_values_ml(ml.net ,"SpType",
                  actors=higher.level[higher.level %in% actors_ml(ml.net)],
                  values="higher.level")
    set_values_ml(ml.net ,"SpType",
                  actors=lower.level[lower.level %in% actors_ml(ml.net)],
                  values="lower.level")
    ## ml.net.all  <- as.igraph(ml.net, merge.actors = FALSE)
    ## ml.net.merged  <- as.igraph(ml.net, merge.actors = TRUE)


    pdf.f(plotMl, sprintf("figures/net%s%s.pdf", yr,
                          paste(species, collapse="")),
          width=10, height=6)


    write_ml(ml.net, file=sprintf("saved/net%s%s.mpx", yr,
                                   paste(species, collapse="")),
             merge.actors = FALSE)

    ## save(ml.net.all, ml.net.merged,
    ##      file=sprintf("saved/multilayer%s%s.Rdata", yr,
    ##                    paste(species, collapse="")))

}


