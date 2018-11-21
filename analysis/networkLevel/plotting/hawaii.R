rm(list=ls())
setwd('~/Dropbox/skyIslands/analysis/networkLevel')
source('src/initialize.R')
load(file='saved/corMets.Rdata')
## ************************************************************


pdf.f <- function(f, file, ...) {
  cat(sprintf('Writing %s\n', file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}


calcCI <- function(dats, colname){
  ci.ub <- dats[, paste("z", colname, sep=".")] +
    dats[, paste("sd", colname, sep=".")]*1.96

  ci.lb <- dats[, paste("z", colname, sep=".")] -
    dats[, paste("sd", colname, sep=".")]*1.96
  return(cbind(ci.ub, ci.lb))
}

plot.panels <- function(){
  f <- function(){
    layout(matrix(1:6, nrow=2, byrow=FALSE))
    par(oma=c(3.5, 3, 0.5, 1),
        mar=c(1, 4, 0.5, 2), cex.axis=1.2)
    plot(cor.dats$zniche.overlap.HL[order(cor.dats$Lat)] ~
             cor.dats$Lat[order(cor.dats$Lat)],
         xlim=range(cor.dats$Lat),
         ylim=range(c(0, cor.dats$zniche.overlap.HL)) + c(0, 0.01),
         ylab="",
         pch=16,
         cex=1.3,
         xlab="",
         xaxt="n")

    mtext("Pollinators", 2,
          line=5, cex=1.5)
    mtext("Niche overlap", 2,
          line=3, cex=1.2)

    plot(cor.dats$zniche.overlap.LL[order(cor.dats$Lat)] ~
             cor.dats$Lat[order(cor.dats$Lat)],
         xlim=range(cor.dats$Lat),
         ylim=range(c(0, cor.dats$zniche.overlap.LL)) + c(0, 0.01),
         ylab="",
         pch=16,
         cex=1.3,
         xlab="")


    mtext("Plants", 2,
          line=5, cex=1.5)
    mtext("Niche overlap", 2,
          line=3.2, cex=1.2)
    mtext("Latitude", 1, line=3, cex=1.2)


    plot(cor.dats$number.of.species.HL[order(cor.dats$Lat)] ~
             cor.dats$Lat[order(cor.dats$Lat)],
         xlim=range(cor.dats$Lat),
         ylim=range(c(0, cor.dats$number.of.species.HL)) + c(0, 10),
         ylab="",
         pch=16,
         cex=1.3,
         xlab="",
         xaxt="n")

    mtext("Richness", 2,
          line=3.2, cex=1.2)

    plot(cor.dats$number.of.species.LL[order(cor.dats$Lat)] ~
             cor.dats$Lat[order(cor.dats$Lat)],
         xlim=range(cor.dats$Lat),
         ylim=range(c(0, cor.dats$number.of.species.LL)) + c(0, 20),
         ylab="",
         pch=16,
         cex=1.3,
         xlab="")

    mtext("Richness", 2,
          line=3.2, cex=1.2)
    mtext("Latitude", 1, line=3, cex=1.2)

    plot(cor.dats$pol.mean.degree[order(cor.dats$Lat)]/
         cor.dats$number.of.species.LL[order(cor.dats$Lat)] ~
             cor.dats$Lat[order(cor.dats$Lat)],
         xlim=range(cor.dats$Lat),
         ylim=range(c(0, cor.dats$pol.mean.degree/cor.dats$number.of.species.LL),
                    na.rm=TRUE),
         ylab="",
         pch=16,
         cex=1.3,
         xlab="",
         xaxt="n")

    mtext("Niche breadth", 2,
          line=3.2, cex=1.2)

    plot(cor.dats$plant.mean.degree[order(cor.dats$Lat)]/
         cor.dats$number.of.species.HL[order(cor.dats$Lat)] ~
             cor.dats$Lat[order(cor.dats$Lat)],
         xlim=range(cor.dats$Lat),
         ylim=range(c(0,
                      cor.dats$plant.mean.degree/cor.dats$number.of.species.HL),
                    na.rm=TRUE),
         ylab="",
         pch=16,
         cex=1.3,
         xlab="")

    mtext("Niche breadth", 2,
          line=3.2, cex=1.2)
    mtext("Latitude", 1, line=3, cex=1.2)



  }
  fig.path <- 'figures'
  pdf.f(f, file=file.path(fig.path,
             sprintf("%s.pdf", "guild")),
        width=8.5, height=5)
}

plot.panels()



plot.panels.network <- function(conservative=TRUE){
  f <- function(){
    layout(matrix(1:2, nrow=2, byrow=FALSE))
    par(oma=c(2, 3, 0.5, 1),
        mar=c(2.5, 4, 0.5, 2), cex.axis=1.2)

    ## ci.nodf <- calcCI(cor.dats, "NODF")
    plot(cor.dats$zNODF[order(cor.dats$Lat)] ~
             cor.dats$Lat[order(cor.dats$Lat)],
         ylab="",
         ylim=range(c(0,cor.dats$zNODF))  + c(-1, 0),
         xlim=range(cor.dats$Lat),
         pch=16,
         xlab="",
         xaxt="n")
    mtext("Nestedness", 2,
          line=3.2, cex=1.2)


    ## arrows(x0=log(cor.dats$Lat),
    ##        y0= ci.nodf[,"ci.lb"],
    ##        y1=ci.nodf[,"ci.ub"],
    ##        angle=90,
    ##        length=0, code=3, col="black", lwd=1.5)

    ## ci.mod <- calcCI(cor.dats, "mod.met.R")
    plot(cor.dats$zmod.met.R[order(cor.dats$Lat)] ~
             cor.dats$Lat[order(cor.dats$Lat)],
         ylab="",
         xlab="",
         xlim=range(cor.dats$Lat),
         ylim=range(c(0, cor.dats$zmod.met.R)) + c(0, 1),
         pch=16)

    ## arrows(x0=log(cor.dats$Lat),
    ##        y0= ci.mod[,"ci.lb"],
    ##        y1=ci.mod[,"ci.ub"],
    ##        angle=90,
    ##        length=0, code=3, col="black", lwd=1.5)

    mtext("Modularity", 2,
          line=3.2, cex=1.2)
    mtext("Latitude", 1, line=3, cex=1.2)
  }
  fig.path <- 'figures'
  pdf.f(f, file=file.path(fig.path,
             sprintf("%s.pdf", "networks")),
        width=4, height=5)
}

plot.panels.network()




