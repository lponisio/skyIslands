##mantel tests
##https://jkzorz.github.io/2019/07/08/mantel-test.html

##Species abundance dissimilarity matrix: created using a distance measure, i.e. Bray-curtis dissimilarity. This is the same type of dissimilarity matrix used when conducting an ANOSIM test or when making an NMDS plot
#Geographic distance matrix: the physical distance between sites (i.e. Haversine distance)

library(vegan)
library(geosphere)
library(betapart)
library(tidyverse)


meta_cols <- c('UniqueID', 'Family', 'Genus', 'Species', 'GenusSpecies', 'Site', 'Lat', 'Long')

spec16s <- spec.net %>%
  filter(Apidae == 1) %>%
  select(all_of(meta_cols), starts_with('16s')) %>%
  na.omit()

##Species abundance dissimilarity matrix: created using a distance measure, i.e. Bray-curtis dissimilarity. This is the same type of dissimilarity matrix used when conducting an ANOSIM test or when making an NMDS plot
#Geographic distance matrix: the physical distance between sites (i.e. Haversine distance)


custom.plot.decay <- function(x, xlim = c(0, max(x$data.x)), ylim = c(0, 1), add = FALSE, 
                              remove.dots = FALSE, col = "black", pch = 1, lty = 1, lwd = 5, 
                              cex = 1, ...) 
{
  if (!inherits(x, "decay")) {
    stop("The input is not a distance-decay model fitted with decay.model().", 
         call. = TRUE)
  }
  data <- data.frame(as.vector(x$data.x), as.vector(x$data.y))
  if (!remove.dots) {
    pch = pch
  }
  else {
    pch = ""
  }
  if (!add) {
    plot(x$data.x, x$data.y, xlim = xlim, ylim = ylim, col = col, pch = pch, cex = cex, 
         ...)
  }
  if (add) {
    points(x$data.x, x$data.y, xlim = xlim, ylim = ylim, col = col, pch = pch, 
           cex = cex, ...)
  }
  switch(x$y.type, similarities = {
    lines(fitted(x$model)[order(x$data.x)] ~ sort(x$data.x), 
          col = col, lty = lty, lwd = lwd, ...)
  }, dissimilarities = {
    lines(1 - fitted(x$model)[order(x$data.x)] ~ sort(x$data.x), 
          col = col, lty = lty, lwd = lwd, ...)
  })
}




## This function takes as input a data frame and a bee genus and filters the data frame to
## just bees of the input genus, computes bray curtis dissimilarity of their 16s communities,
## computes haversine distance matrix of sample sites, performs a mantel test to determine
## correlation between community dissimilarity and distance, then plots the distance decay curves

genusspecies.decay.model <- function(data, type, which){
  #bray curtis dissimilarity matrix of 16s
  
  
  if(type == 'Genus'){
    abund <- data %>%
      filter(Genus == which) %>%
      select(UniqueID, starts_with('16s')) %>%
      select(-UniqueID)
    
    #distance matrix of sites
    geo <- data %>%
      filter(Genus == which) %>%
      select(UniqueID, Long, Lat) %>%
      select(-UniqueID) %>%
      mutate()
  } else if (type == 'GenusSpecies'){
    
    abund <- data %>%
      filter(GenusSpecies == which) %>%
      select(UniqueID, starts_with('16s')) %>%
      select(-UniqueID)
    #distance matrix of sites
    geo <- data %>%
      filter(GenusSpecies == which) %>%
      select(UniqueID, Long, Lat) %>%
      select(-UniqueID) %>%
      mutate()
  }
  
  
  
  
  #abundance data frame - bray curtis dissimilarity
  dist.abund <- vegdist(abund, method = "bray")
  
  #geographic data frame - haversine distance in m (takes a df with lat and long and calculates dist)
  d.geo <- distm(geo, fun = distHaversine)
  dist.geo <- as.dist(d.geo)/1000
  
  #abundance vs geographic mantel test
  abund_geo  = mantel(dist.abund, dist.geo, method = "spearman", permutations = 9999, na.rm = TRUE)
  print(abund_geo)
  
  dist_decay_model <- betapart::decay.model(dist.abund,
                                            dist.geo,
                                            y.type='dissim',
                                            model.type = 'exp',
                                            perm=100)
  # dist_decay_plot <- plot.decay(dist_decay_model,
  #                               main=genus)
  # dist_decay_plot
  dist_decay_model
  
}


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
