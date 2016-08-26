##the purpose of this function is to break up data with many
##sites/years and prepare it for network analysis.

break.net <- function(spec.dat,
                      by.site = TRUE,
                      by.year = FALSE,
                      plant.nam = "PlantGenSp",
                      pol.nam = "GenSp"){
  samp2site.spp <- function(site, spp, abund) { 
    x <- tapply(abund, list(site= site, spp= spp), sum)
    x[is.na(x)] <- 0
    return(x)
  }
  drop.net <- function(z){
    z[!sapply(z, FUN=function(q){
      any(dim(q) < 3)
    })]
  }

  ## status <- split(spec.dat, spec.dat$SiteStatus)
  if(by.site){
    networks <- split(spec.dat, spec.dat$Site)
  }
  if(by.year){
    networks <- lapply(networks, function(x){
      lapply(split(x, f=x$Year), as.matrix)
    })
  }
  ## formats data matrices appropriate for network analysis
  if(by.site & by.year){
    comms <- rapply(networks, function(y){
      samp2site.spp(site=y[, plant.nam],
                    spp=y[, pol.nam],
                    abund=rep(1, nrow(y)))
    }, how="replace")
    adj.mat <- unlist(lapply(comms, drop.net), recursive=FALSE)
  }else{
    comms <- lapply(networks, function(y){
      samp2site.spp(site=y[, plant.nam],
                    spp=y[, pol.nam],
                    abund=rep(1, nrow(y)))
    })
    adj.mat <- comms
  }
  return(adj.mat) 
}
