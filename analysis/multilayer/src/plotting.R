
plotMl <- function(){
    sp.types <- get_values_ml(ml.net,
                              actors = vertices_ml(ml.net)[[1]],
                              attribute = "SpType")
    gr <- values2graphics(sp.types[[1]])

    ## ed <- edges_ml(ml.net)
    ## ed.type <- rep("intralayer", nrow(ed))
    ## ed.type[ed$from_layer != ed$to_layer] <- "interlayer"

    ## gr.ed <- values2graphics(ed.type)

    l <- layout_multiforce_ml(ml.net, w_inter = 1,
                              gravity = 1)
    bk <- par("mar")
    par(mar=c(0,0,0,0))
    plot(ml.net, layers=spatial.order[[yr]],
         layout = l, # vertex.labels = "",
         vertex.labels.cex =0.5,
         edge.col="grey",
         vertex.color = gr$color)
    par(mar=bk)

    legend("bottomright",
           legend = gr$legend.text,
           col = gr$legend.col,
           pt.bg = gr$legend.col,
           pch = gr$legend.pch,
           bty = "n", pt.cex = 1, cex = .5,
           inset = c(0.05, 0.05)
           )
}
