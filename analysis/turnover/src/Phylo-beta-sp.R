

f.curve <- function(x, intercept, slope){
  y <- intercept + x * slope
  return(y)
}

## phylo beta of plants, pol

pp.dis <- function(spec.dat, types, c.dist, sub, path){
  comm.dis <- function(spec.dat, type){
    prep.comm <- aggregate(spec.dat[, type[1]],
                           list(sp=spec.dat[, type[1]],
                                sp2= spec.dat[, type[2]]), length)
    comm <- samp2site.spp(prep.comm$sp, prep.comm$sp2, prep.comm$x)
    comm.dis <- vegdist(comm, "chao", diag= TRUE)
    dengram <- hclust(comm.dis, method="average")
    fden <- function(){
      plot(dengram, cex.lab=0.5)
    }
    pdf.f(fden, file= file.path(path, sprintf("%s.pdf",
                  paste("dedrogram", type[1], sep="_"))),
          width=25, height=5)
    
    prep.site.comm <- aggregate(spec.dat[, type[1]],
                                list(site=spec.dat$Site,
                                     sp= spec.dat[, type[1]]), length)
    site.comm <- samp2site.spp(prep.site.comm$site,
                               prep.site.comm$sp,
                               prep.site.comm$x)
    phylo.b <- as.matrix(comdist(site.comm, cophenetic(dengram),
                                 abundance.weighted=TRUE))
    return(phylo.b[lower.tri(phylo.b)])
  }
  ##pollinator dissimilarity
  c.pol <-comm.dis(spec.dat, types[1:2])
  ##plant dissimilarity
  c.plant <- comm.dis(spec.dat, types[2:1])
  c.dist <- scale(log(c.dist))
  ##linear models
  lm.pol <- lm(c.pol~c.dist)
  lm.plant <- lm(c.plant~c.dist)

  f <- function(){
    plot(c.pol~c.dist, pch =16, col="darkolivegreen", ylim=c(0,1),
         xlab="distance km (scaled)", ylab="dissimilarity")
    points(c.plant~c.dist, pch=16, col="dodgerblue")
    legend("bottomright", legend= c("Pollinators", "Plants"),
           col=c("darkolivegreen", "dodgerblue"), pch=16, bty="n", cex=1.5)
    curve(f.curve(x=x, intercept =
                  coef(summary(lm.pol))[1,1],
                  slope=coef(summary(lm.pol))[2,1]),
          lwd=2, col="darkolivegreen", add=TRUE)
    curve(f.curve(x=x, intercept =
                  coef(summary(lm.plant))[1,1],
                  slope=coef(summary(lm.plant))[2,1]),
          lwd=2, col="dodgerblue", add=TRUE)
  } 
  pdf.f(f, file= file.path(path, sprintf("%s.pdf",
             paste("phylo-beta", types[1], sub, sep="_"))),
        width=5, height=5)
  
  return(list(summary(lm.pol), summary(lm.plant))) 
}
