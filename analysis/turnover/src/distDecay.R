
ppint.dis <- function(spec.dat,
                      types= c("GenusSpecies", "PlantGenusSpecies", "Int"),
                      c.dist, sub="all", method.dist="gower"){
    comm.dis <- function(spec.dat, type){
        browser()
        prep.comm <- aggregate(spec.dat[, type],
                               list(Site=spec.dat$Site,
                                    Year=spec.dat$Year,
                                    GenusSpecies=spec.dat[, type]),
                               length)
        prep.comm$Nyears  <- 20
        prep.comm$Nyears[prep.comm$Site %in% c("HM", "SM", "PL", "VC")] <- 10
        prep.comm$AveAbund <- prep.comm$x/prep.comm$Nyears
        comm <- samp2site.spp(prep.comm$Site,
                              prep.comm$GenusSpecies, prep.comm$AveAbund)
        comm.dis <- as.matrix(vegdist(comm, method.dist, diag= TRUE))
        c.comm <- comm.dis[lower.tri(comm.dis)]
        return(c.comm)
    }
    browser()
    ##pollinator dissimilarity
    c.pol <-comm.dis(spec.dat, types[1])
    ##plant dissimilarity
    c.plant <- comm.dis(spec.dat, types[2])
    ##interaction dissimilarity
    c.int <-  comm.dis(spec.dat, types[3])
    c.dist <- scale(log(c.dist))
    ##linear models
    lm.pol <- lm(c.pol~c.dist)
    lm.plant <- lm(c.plant~c.dist)
    lm.int <- lm(c.int~c.dist)

    return(list(summary(lm.pol), summary(lm.plant), summary(lm.int)))
}




f <- function(){

    f.curve <- function(x, intercept, slope){
        y <- intercept + x * slope
        return(y)
    }
    layout(1)
    plot(c.pol~c.dist, pch =16, col="darkolivegreen", ylim=c(0,1),
         xlab="distance", ylab="dissimilarity")
    points(c.plant~c.dist, pch=16, col="dodgerblue")
    points(c.int~c.dist, pch=16, col="darkred")
    legend("bottomright", legend= c("Pollinators", "Plants",
                                    "Interactions"), col=c("darkolivegreen",
                                                           "dodgerblue","darkred"),
           pch=16, bty="n")

    curve(f.curve(x=x, intercept =
                           coef(summary(lm.pol))[1,1],
                  slope=coef(summary(lm.pol))[2,1]),
          lwd=2, col="darkolivegreen", add=TRUE)
    curve(f.curve(x=x, intercept =
                           coef(summary(lm.plant))[1,1],
                  slope=coef(summary(lm.plant))[2,1]),
          lwd=2, col="dodgerblue", add=TRUE)

    curve(f.curve(x=x, intercept =
                           coef(summary(lm.int))[1,1],
                  slope=coef(summary(lm.int))[2,1]),
          lwd=2, col="darkred", add=TRUE)

}
path <- '~/Dropbox/SkyIslands/analysis/figures/distanceDecay'
pdf.f(f, file= file.path(path, sprintf("%s.pdf",
                                       paste(method.dist, "dd", types[1], sub, sep="_"))),
      width=5, height=5)
