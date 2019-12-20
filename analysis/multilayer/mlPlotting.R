library(multinet)

net.2012 <- read_ml(file="saved/net2012.mpx")

## ********************************************************************
## plotting
## ********************************************************************
sp.types <- get_values_ml(net.2012,  actors = vertices_ml(net.2012)[[1]],
                          attribute = "SpType")
gr <- values2graphics(sp.types[[1]])

## ********************************************************************
## modularity
## ********************************************************************

ml_clust <- glouvain_ml(net.2012)
ml_abacus <- abacus_ml(net.2012, min.actors=3, min.layers=1)
ml_cp <- clique_percolation_ml(net.2012, k=2, m=1)
ml_infomap <- infomap_ml(net.2012)


l <- layout_multiforce_ml(net.2012)

f1 <- function(){
    plot(net.2012,
         com = ml_clust,
         layers=c("JC", "SC", "MM", "PL",  "CH"),
         vertex.labels = "",
         layout=l, grid = c(1,5),
         vertex.color = gr$color
         )
}

pdf.f(f1, "figures/glouvain.pdf", width=10, height=6)

f2 <- function(){
    plot(net.2012,
         com = ml_abacus,
         layers=c("JC", "SC", "MM", "PL",  "CH"),
         vertex.labels = "",
         layout=l, grid = c(1,5),
         vertex.color = gr$color
         )
}

pdf.f(f2, "figures/abacus.pdf", width=10, height=6)

f3 <- function(){
    plot(net.2012,,
         com = ml_cp,
         layers=c("JC", "SC", "MM", "PL",  "CH"),
         vertex.labels = "",
         layout=l, grid = c(1,5),
         vertex.color = gr$color
         )
}
pdf.f(f3, "figures/cp.pdf", width=10, height=6)


f4 <- function(){
    plot(net.2012,
         com = ml_infomap,
         layers=c("JC", "SC", "MM", "PL",  "CH"),
         vertex.labels = "",
         layout=l, grid = c(1,5),
         vertex.color = gr$color
         )
}
pdf.f(f4, "figures/infomap.pdf", width=10, height=6)


