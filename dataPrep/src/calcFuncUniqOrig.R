## Functional roles (Coux et al. 2016) and community functional diversity.
## All community scales share one species pool, Gower matrix and PCoA.

## A community matrix with explicit rows for sampled communities, including
## those without bees or without taxa represented in the trait data.
make.comm.matrix <- function(spec, events, keys, taxa, suffix) {
  groups <- events %>% distinct(across(all_of(keys))) %>%
    arrange(across(all_of(keys))) %>%
    mutate(commID=paste0(suffix, '.', row_number()))
  counts <- spec %>% select(all_of(keys), GenusSpecies) %>%
    left_join(groups, by=keys, relationship='many-to-one')
  a <- with(counts, table(factor(commID, levels=groups$commID),
                          factor(GenusSpecies, levels=taxa)))
  a <- matrix(as.numeric(a), nrow=nrow(groups), ncol=length(taxa),
              dimnames=list(groups$commID, taxa))
  coverage <- counts %>% group_by(commID) %>% summarise(
    BeeAbundance=n(),
    BeeIdentifiedAbundance=sum(!is.na(GenusSpecies) & GenusSpecies != ''),
    BeeIdentifiedRichness=n_distinct(GenusSpecies[!is.na(GenusSpecies) &
                                                   GenusSpecies != '']),
    .groups='drop')
  groups <- groups %>% left_join(coverage, by='commID') %>%
    mutate(across(c(BeeAbundance, BeeIdentifiedAbundance,
                    BeeIdentifiedRichness), ~replace_na(.x, 0)),
           BeeTraitAbundance=as.numeric(rowSums(a)),
           BeeTraitRichness=as.numeric(rowSums(a > 0)),
           BeeTraitCoverage=if_else(BeeAbundance > 0,
                                    BeeTraitAbundance/BeeAbundance, NA_real_))
  list(groups=groups, a=a)
}

## Species-level raw distances, not within-community z-scores.
## Uniqueness is independent of abundance; originality uses the weighted
## centroid, as Coux et al. did for pollinators (w.abun=TRUE).
## Also retain the unweighted species-centered version.
calc.species.roles <- function(coords, abundance, w.abun=TRUE) {
  n <- nrow(coords)
  p <- if (w.abun) abundance/sum(abundance) else rep(1/n, n)
  center <- colSums(coords*p)
  originality <- sqrt(rowSums(sweep(coords, 2, center)^2))
  unweightedOriginality <- sqrt(rowSums(sweep(coords, 2, colMeans(coords))^2))
  uniq <- rep(NA_real_, n)
  if (n > 1) {
    d <- as.matrix(dist(coords))
    diag(d) <- Inf
    uniq <- apply(d, 1, min)
  }
  data.frame(GenusSpecies=rownames(coords), uniq=unname(uniq),
             originality=unname(originality),
             unweightedOriginality=unname(unweightedOriginality), row.names=NULL)
}

calcFuncUniqOrigComm <- function(traits, weights, a, w.abun=TRUE,
                                corr='cailliez') {
  a <- as.matrix(a[, rownames(traits), drop=FALSE])
  fd <- data.frame(commID=rownames(a), BeeFDis=NA_real_,
                   BeeFEve=NA_real_, BeeFDiv=NA_real_)
  empty.roles <- data.frame(commID=character(), GenusSpecies=character(),
                            uniq=numeric(), originality=numeric(),
                            unweightedOriginality=numeric())
  keep <- rowSums(a) > 0
  if (!any(keep)) return(list(fd=fd, by.comm.mets=empty.roles,
                              coords=NULL, distances=NULL, eigenvalues=NULL))
  present.a <- a[keep, , drop=FALSE]
  quality <- NA_real_
  if (nrow(traits) == 1) {
    trait.dist <- dist(matrix(0, 1, 1, dimnames=list(rownames(traits), NULL)))
    coords <- matrix(0, 1, 1, dimnames=list(rownames(traits), 'Axis1'))
    eigenvalues <- 0
    fd$BeeFDis[keep] <- 0
  } else {
    ## Gower range-standardizes ITD and treats categorical traits as nominal.
    ## Missing traits are omitted pairwise, not coded as a new category.
    trait.dist <- FD::gowdis(traits, w=weights)
    if (anyNA(trait.dist)) {
      bad <- which(is.na(as.matrix(trait.dist)), arr.ind=TRUE)
      stop('No comparable traits for some species pairs: ',
           paste(head(rownames(traits)[bad[,1]], 10), collapse=', '),
           '. Review the trait data rather than assigning arbitrary distances.')
    }
    if (all(trait.dist == 0)) {
      coords <- matrix(0, nrow(traits), 1,
                       dimnames=list(rownames(traits), 'Axis1'))
      eigenvalues <- 0
      fd$BeeFDis[keep] <- 0
    } else {
      ## dbFD uses a single regional ordination. Its default m="max"
      ## selects the number of axes for FRic and FDiv from community richness.
      ## FRic must be enabled to calculate FDiv, but need not be exported.
      calc.hull <- any(rowSums(present.a > 0) >= 3)
      fit <- FD::dbFD(trait.dist, a=present.a, w.abun=w.abun,
                      corr=corr, calc.FRic=calc.hull, calc.FDiv=calc.hull,
                      calc.CWM=FALSE, print.pco=TRUE)
      coords <- as.matrix(fit$x.axes)
      eigenvalues <- fit$x.values
      if (!is.null(fit$qual.FRic)) quality <- fit$qual.FRic
      fd$BeeFDis[keep] <- unname(fit$FDis)
      fd$BeeFEve[keep] <- unname(fit$FEve)
      if (!is.null(fit$FDiv)) fd$BeeFDiv[keep] <- unname(fit$FDiv)
    }
  }
  roles <- lapply(which(keep), function(i) {
    present <- a[i, ] > 0
    x <- coords[colnames(a)[present], , drop=FALSE]
    abundance <- as.numeric(a[i, present])
    out <- calc.species.roles(x, abundance, w.abun=w.abun)
    out$commID <- rownames(a)[i]
    out
  })
  fd <- fd %>% mutate(across(c(BeeFDis, BeeFEve, BeeFDiv),
                             ~replace(.x, !is.finite(.x), NA_real_)))
  list(fd=fd, by.comm.mets=bind_rows(roles), coords=coords,
       distances=trait.dist, eigenvalues=eigenvalues, qual.FRic=quality)
}
