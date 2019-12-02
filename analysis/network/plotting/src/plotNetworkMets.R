
plotNetworkMets <- function(){
    col.lines <- rep("black", length(ys))
    col.fill <- add.alpha(col.lines, alpha=0.3)
    names(col.lines) <- names(col.fill) <- ys


    layout(matrix(1:length(ys), ncol=2, byrow=TRUE))


    par(oma=c(6, 3, 2, 1),
        mar=c(0.5, 7, 1, 1), cex.axis=1.5)


    for(y in ys){
        print(y)
        ## create a matrix of possible variable values
        dd.met <- expand.grid(xvar =seq(
                                  from=  min(cor.dats[, x]),
                                  to= max(cor.dats[, x]),
                                  length=10),
                              Year=c("2012", "2017", "2018"),
                              yvar=0)
        colnames(dd.met)[1] <- x
        colnames(dd.met)[3] <- y

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
                   plot.x=FALSE,
                   years=unique(cor.dats$Year))
        if(grepl("number.of.species", y)){
            axis(1, pretty(cor.dats[,x], 4))
            mtext(xlabel, 1, line=3, cex=1)
        }

        plotDiag <- function(){
            plotDiagnostics(mods=mods.div[[y]], dats=cor.dats)
        }

        pdf.f(plotDiag,
              file=file.path('figures/diagnostics/',
                             sprintf('%s_%s.pdf', y, xlabel)),
              height=6, width=4)
    }

}

