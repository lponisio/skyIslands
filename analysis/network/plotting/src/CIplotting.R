plot.panel <- function(dats,
                       new.dd,
                       y1,
                       y2,
                       xs,
                       col.lines,
                       col.fill,
                       ylabel,
                       plot.x=TRUE,
                       FUN=mean,
                       pchs=c(16),
                       years){
    plotting.loop <- function(){
        for(yr in years){
            ## take means of ys for each plot
            ## ys <- aggregate(list(y=dats[,y1]),
            ##                 list(x=dats[,xs]),
            ##                 mean, na.rm=TRUE)

            these.dats <- dats[dats$Year == yr,]
            pred.dats <- new.dd[new.dd$Year == yr,]
            ys <- data.frame(y=these.dats[,y1], x=these.dats[,xs])

            ## plots means
            points(x=ys$x,
                   y=ys$y,
                   pch=pchs[yr],
                   col=col.lines[y1],
                   cex=1.2)

            if(nrow(pred.dats) > 1){
                ## add fill from ci.up to ci.lb
                polygon(c(pred.dats[,xs],
                          rev(pred.dats[,xs])),
                        c(pred.dats$phi,
                          rev(pred.dats$plo)),
                        col=col.fill[y1], border=NA)
                ## plots CI
                lines(x=pred.dats[,xs],
                      y=pred.dats[,y1],
                      col=col.lines[y1],
                      lwd=2)
                lines(x=pred.dats[,xs],
                      y=pred.dats$plo,
                      col=col.lines[y1],
                      lty="dashed")
                lines(x=pred.dats[,xs],
                      y=pred.dats$phi,
                      col=col.lines[y1],
                      lty="dashed")
            }
        }
    }
    plot(NA, xlim=range(dats[, xs], na.rm=TRUE),
         ylim=range(c(new.dd$plo, new.dd$phi,
                      dats[, y1]) ,
                    na.rm=TRUE),
         xlab="",
         ylab="",
         xaxt="n",
         yaxt="n",
         las=1)

    axis(2, pretty(range(c(0, new.dd$plo, new.dd$phi,
                           dats[, y1]),
                         na.rm=TRUE),
                   n = 5),
         las=1)
    mtext(ylabel, 2, line=4, cex=1)

    if(plot.x){
        axis(1, pretty(dats[,xs], 4))
    }
    plotting.loop()
}



plot.panel_nofill <- function(dats,
                       new.dd,
                       y1,
                       y2,
                       xs,
                       col.lines,
                       pchs=c(16),
                       years){
    plotting.loop <- function(){
        for(yr in years){
            ## take means of ys for each plot

            these.dats <- dats[dats$Year == yr,]
            pred.dats <- new.dd[new.dd$Year == yr,]
            ys <- data.frame(y=these.dats[,y1], x=these.dats[,xs])

            ## plots means
            points(x=ys$x,
                   y=ys$y,
                   pch=pchs[yr],
                   col="black",
                   cex=1.2)

            if(nrow(pred.dats) > 1){
                ## add fill from ci.up to ci.lb
                ## plots CI
                lines(x=pred.dats[,xs],
                      y=pred.dats[,y1],
                      col=col.lines[y1],
                      lwd=2)
            }
        }
    }
    plotting.loop()
}
