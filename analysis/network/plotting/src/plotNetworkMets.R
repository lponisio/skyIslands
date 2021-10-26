library(viridis)
plotNetworkMetsPlantPol <- function(){
    col.lines <- rep("black", length(ys))
    col.fill <- add.alpha(col.lines, alpha=0.3)
    names(col.lines) <- names(col.fill) <- ys
    col.fill[c("niche.overlap.LL",
               "zweighted.NODF",
               "number.of.species.LL",
               "zH2", "mean.number.of.links.LL",
               "mean.number.of.links.HL")] <-
        add.alpha("white", alpha=0.1)

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
                              Doy=mean(cor.dats$Doy),
                              yvar=0)
        colnames(dd.met)[1] <- x
        colnames(dd.met)[4] <- y

        if(species[2] == "Parasite"){
            dd.met$Year <- NULL
            }


        ## "predict" data using model coefficient to draw predicted
        ## relationship and CI
        met.pi <- predict.int(mod= mods.div[[y]],
                              dd=dd.met,
                              y=y,
                              family="guassian")

        if(species[1] == "Plant" & length(years[[y]]) == 1){
            print(paste("one year", y))
            met.pi <- met.pi[met.pi$Year == yrs[3],]
        }

        pchs <- c(16, 15, 17)
        names(pchs) <- yrs

        cor.dats$SiteStatus <- "all"
        plot.panel(dats=cor.dats,
                   new.dd=met.pi,
                   xs=x,
                   y1=y,
                   col.lines=col.lines,
                   col.fill=col.fill,
                   ylabel= ylabs[y],
                   plot.x=FALSE,
                   years=yrs,
                   pch=pchs)

        if(grepl("NODF", y)|grepl("H2", y)){
            axis(1, pretty(cor.dats[,x], 4))
            mtext(xlabel, 1, line=3, cex=1)
        }

        if(grepl("niche.overlap.HL", y)){
            legend("topright", legend=c("2012", "2017", "2018"),
                   pch=pchs, bty="n", cex=1.2)
        }

        plotDiag <- function(){
            plotDiagnostics(mods=mods.div[[y]], dats=cor.dats)
        }

        pdf.f(plotDiag,
              file=file.path('figures/diagnostics',
                             sprintf('%s_%s.pdf', y, xlabel)),
              height=6, width=4)
    }

}



plotNetworkMetsPlantPolDOY <- function(){
    x <- "Doy"

    layout(matrix(1:length(ys), ncol=2, byrow=TRUE))
    par(oma=c(6, 3, 2, 1),
        mar=c(0.5, 7, 1, 1), cex.axis=1.5)


    yrs <- c("2012", "2017", "2018")
    years <- list(niche.overlap.LL=yrs[1],
                  niche.overlap.HL=yrs[1],
                  weighted.cluster.coefficient.LL=yrs,
                  weighted.cluster.coefficient.HL=yrs,
                  mean.number.of.links.LL=yrs,
                  mean.number.of.links.HL=yrs,
                  number.of.species.LL=yrs,
                  number.of.species.HL=yrs,
                  zweighted.NODF=yrs,
                  H2=yrs)

    cols <- magma(length(unique(cor.dats$Lat)))
    for(y in ys){
        print(y)

        ylims <- c(min(cor.dats[, y], na.rm=TRUE) -
                   abs(mean(cor.dats[, y], na.rm=TRUE)),
                   max(cor.dats[, y], na.rm=TRUE) +
                   abs(mean(cor.dats[, y], na.rm=TRUE)))

        plot(NA, xlim=range(cor.dats[, x], na.rm=TRUE),
             ylim=ylims,
             xlab="",
             ylab="",
             xaxt="n",
             yaxt="n",
             las=1)

        axis(2, pretty(ylims,
                       n = 5),
             las=1)
        mtext(ylabs[y], 2, line=4, cex=1)
        if(grepl("niche.overlap.HL", y)){
            legend("topright", legend=c("2012", "2017", "2018"),
                   pch=pchs, bty="n", cex=1.2)
        }


        if(grepl("NODF", y)|grepl("H2", y)){
            axis(1, pretty(cor.dats[,x], 4))
            mtext("Day of the year", 1, line=3, cex=1)
        }

        for(j in 1:length(unique(cor.dats$Lat))){
            ## create a matrix of possible variable values
            dd.met.mean <- expand.grid(Doy =seq(
                                           from=  min(cor.dats$Doy),
                                           to= max(cor.dats$Doy),
                                           length=10),
                                       Year=c("2012", "2017", "2018"),
                                       Lat=unique(cor.dats$Lat)[j],
                                       yvar=0)
            colnames(dd.met.mean)[4] <- y

            if(species[2] == "Parasite"){
                dd.met.mean$Year <- "2018"
            }

            ## "predict" data using model coefficient to draw predicted
            ## relationship and CI
            met.pi.mean <- predict.int(mod= mods.div[[y]],
                                       dd=dd.met.mean,
                                       y=y,
                                       family="guassian")

            if(length(years[[y]]) == 1){
                print(paste("one year", y))
                met.pi.mean <- met.pi.mean[met.pi.mean$Year == yrs[3],]

            }

            pchs <- c(16, 15, 17)
            names(pchs) <- yrs

            col.lines.mean <- rep(cols[j], length(ys))

            names(col.lines.mean) <- ys

            plot.panel_nofill(dats=cor.dats,
                              new.dd=met.pi.mean,
                              xs=x,
                              y1=y,
                              col.lines=col.lines.mean,
                              years=yrs,
                              pch=pchs)

        }
    }

}


