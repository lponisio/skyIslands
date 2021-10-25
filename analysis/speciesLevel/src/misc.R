

pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}

samp2site.spp <- function(site,spp,abund) {
  x <- tapply(abund, list(site=site,spp=spp), sum)
  x[is.na(x)] <- 0

  return(x)
}
