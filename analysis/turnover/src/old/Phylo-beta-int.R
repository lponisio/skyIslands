##phylo-beta of interactions

int.dis <- function(spec.dat, c.dist, abund.w, path){
  comm.dis <- function(spec.dat, abund.w){
    load('src/lc.Rdata')
    prep.site.comm.int <- aggregate(spec.dat[, "Int"],
                                    list(site=spec.dat$Site,
                                         sp= spec.dat[, "Int"]), length)
    site.comm.int <- samp2site.spp(prep.site.comm.int$site,
                                   prep.site.comm.int$sp,
                                   prep.site.comm.int$x)
    node.num <- cbind(lc$hclust$labels, paste(lc$edgelist[,1],
                                              lc$edgelist[,2]))
    node.order <- node.num[,1][match(colnames(site.comm.int),
                                     node.num[,2])]
    colnames(site.comm.int) <- node.order
    phylo.b <- comdist(site.comm.int, cophenetic(lc$hclust),
                       abundance.weighted=abund.w)
    return(phylo.b[lower.tri(phylo.b)])
  }
  c.int <-comm.dis(spec.dat, abund.w)

  ##linear models
  c.dist <- scale(log(c.dist))
  lm.int <- lm(c.int ~ c.dist)
  
  f.curve <- function(x, intercept, slope){
    y <- intercept + x * slope
    return(y)
  }
  f <- function(){
    layout(1)
    plot(c.int~c.dist, pch =16, col="darkred", ylim=c(0,1),
         xlab="distance km (scaled)", ylab="dissimilarity")
    curve(f.curve(x=x, intercept =
                  coef(summary(lm.int))[1,1],
                  slope=coef(summary(lm.int))[2,1]),
          lwd=2, col="darkred", add=TRUE)
  }
  pdf.f(f, file= file.path(path, sprintf("%s.pdf",
             paste("phylo-beta", "int", sep="_"))),
        width=5, height=5)
  return(summary(lm.int))
}
