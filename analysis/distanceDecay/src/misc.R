mut.adj <- function(x) {
  nr <- dim(x)[1]
  nc <- dim(x)[2]
  to.fill <- matrix(0, ncol=nc + nr, nrow=nc + nr)
  to.fill[1:nr,(nr+1):(nc+nr)] <- x
  rownames(to.fill) <- c(rownames(x), rep(NA, nrow(to.fill)-nrow(x)))
  colnames(to.fill) <- c(colnames(x), rep(NA, ncol(to.fill)-ncol(x)))
  adj.mat <- graph.adjacency(to.fill, mode= "undirected",
                             weighted=TRUE,
                             add.rownames="code", add.colnames="names")
  return(adj.mat)
}

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

## This functions takes site-species-abundance data and creates a
## matrix where the sites are columns and the rows are species. 

##this function takes three vectors: 1) the site species are found in.
## For plant animal networks this is generally the pollinators
## 2) the species found in those sites.
## In the case of p-a data is this usally the plants and
## 3) the adbundance of those species. In p-a data this is the
## frequency of interactions 

