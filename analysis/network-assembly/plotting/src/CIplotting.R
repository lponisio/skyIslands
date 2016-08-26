plot.panel <- function(dats,
                       new.dd,
                       y1,
                       xs,
                       col.lines,
                       col.fill,
                       plotPoints=TRUE,...){
  plotting.loop <- function(){
    ys <- aggregate(list(y=dats[,y1]),
                    list(Site=dats$Site,
                         x=dats[,xs]),
                    mean, na.rm=TRUE)
    if(plotPoints){
      points(x=jitter(ys$x, factor=0.25),
             y=ys$y,
             pch=16,
             col=col.fill,
             cex=1.2)
    }
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
    polygon(c(new.dd[,xs],
              rev(new.dd[,xs])),
            c(new.dd$phi,
              rev(new.dd$plo)),
            col=col.fill, border=NA)
  }
  plot(NA,
       xlim=range(new.dd[,xs], na.rm=TRUE),
       ylim=range(c(0, new.dd$phi,  new.dd$plo,
         na.rm=TRUE)),
       xlab="",
       ylab="",
       xaxt="n",
       las=1,...)
  plotting.loop()
}


plot.predict.ypr <- function(new.dd,
                             ylabel,
                             dats,
                             y1,
                             xs="ypr",
                             plotPoints=TRUE,
                             xlabel= "Years Post Restoration"){
  plot.ci <- function(){
    col.lines <-  brewer.pal(3, "Dark2")[2]
    col.fill <- add.alpha(col.lines, alpha=0.2)
    layout(matrix(1, ncol=1))
    par(oma=c(6, 5, 2, 1),
        mar=c(0.5, 0, 0.5, 1))
    plot.panel(dats, new.dd, y1, xs,
               col.lines, col.fill, plotPoints)
    axis(1, pretty(dats[,xs]), labels=pretty(dats[,xs]))
    mtext(ylabel, 2, line=3.5, cex=1.5)
    mtext(xlabel, 1, line=3.5, cex=1.5)
  }
  path <- 'figures' 
  pdf.f(plot.ci, file=file.path(path,
                   sprintf("%s.pdf", paste(
                     gsub(" ", "", ylabel),
                     xlabel, sep="_"))),
        width=4, height=4)

}
