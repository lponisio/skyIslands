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
                       pchs=c(16)){
    plotting.loop <- function(){
        ## take means of ys for each plot
        ys <- aggregate(list(y=dats[,y1]),
                        list(x=dats[,xs]),
                        mean, na.rm=TRUE)

        ## ys <- data.frame(y=dats[,y1], x=dats[,xs])

        ## add fill from ci.up to ci.lb
        polygon(c(new.dd[,xs],
                  rev(new.dd[,xs])),
                c(new.dd$phi,
                  rev(new.dd$plo)),
                col=col.fill, border=NA)
        ## plots means
        points(x=jitter(ys$x, factor=0.25),
               y=ys$y,
               pch=pchs,
               col=col.lines,
               cex=1)
        ## plots CI
        lines(x=new.dd[,xs],
              y=new.dd[,y1],
              col=col.lines,
              lwd=2)
        lines(x=new.dd[,xs],
              y=new.dd$plo,
              col=col.lines,
              lty="dashed")
        lines(x=new.dd[,xs],
              y=new.dd$phi,
              col=col.lines,
              lty="dashed")
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

    axis(2, pretty(range(c(0, 1, new.dd$plo, new.dd$phi,
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
