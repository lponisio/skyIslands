
plotNetworkMets <- function(){
    col.lines <- rep("black", length(ys))
    col.fill <- add.alpha(col.lines, alpha=0.3)
    names(col.lines) <- names(col.fill) <- ys

    layout(matrix(1:length(ys), ncol=2, byrow=TRUE))
    par(oma=c(6, 7, 2, 1),
        mar=c(0.5, 0, 1, 1), cex.axis=1.5)

    for(y in ys){
        print(y)
        browser()
        ## create a matrix of possible variable values
        dd.met <- expand.grid(xvar =seq(
                                  from=  min(scale(cor.dats[, x])),
                                  to= max(scale(cor.dats[, x])),
                                  length=10),
                              yvar=0)
        colnames(dd.met)[1] <- x
        colnames(dd.met)[2] <- y

        ## "predict" data using model coefficient to draw predicted
        ## relationship and CI
        met.pi <- predict.int(mod= mods.div[[y]],
                              dd=dd.met,
                              y=y,
                              family="guassian")

        cor.dats$SiteStatus <- "all"
        plot.panel(dats=cor.dats,
                   new.dd=met.pi,
                   xs=x,
                   y1=y,
                   col.lines=col.lines,
                   col.fill=col.fill,
                   ylabel= ylabs[y],
                   plot.x=FALSE)
        ## if(y == "links.per.species"){
        ##     axis(1, pretty(cor.dats[,x], 4))
        ##     mtext(xlabel, 1, line=3, cex=1)
        ## }

        plotDiag <- function(){
            plotDiagnostics(mods=mods.div[[y]], dats=cor.dats)
        }

        pdf.f(plotDiag,
              file=file.path('figures/diagnostics/',
                             sprintf('%s.pdf', y)),
              height=6, width=4)
    }

}

