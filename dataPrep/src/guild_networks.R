## Per-event plant-pollinator interaction matrices, with explicit guilds and
## transparent identification audits. No annual means and no averaging over
## repeated sampling rounds. Source only after src/specialization.R.

pollinator_guilds <- function(bee.families, butterfly.hawk.families) {
  bs <- bee.families
  lh <- butterfly.hawk.families
  combined <- c(bs, 'Syrphidae', 'Bombyliidae', lh)
  list(
    AllNetted = function(x) rep(TRUE, nrow(x)),
    AllPollinatorGuilds = function(x) x$Family %in% combined,
    Bees = function(x) x$Family %in% bs,
    HoneyBees = function(x) x$GenusSpecies == 'Apis mellifera' &
      x$Family %in% bs,
    BumbleBees = function(x) x$Genus == 'Bombus' & x$Family %in% bs,
    NonBombusBees = function(x) x$Family %in% bs &
      x$Genus != 'Bombus' & x$Genus != 'Apis',
    Syrphids = function(x) x$Family == 'Syrphidae',
    BeeFlies = function(x) x$Family == 'Bombyliidae',
    ButterfliesHawkMoths = function(x) x$Family %in% lh,
    BeesSyrphids = function(x) x$Family %in% c(bs, 'Syrphidae'),
    BeesSyrphidsButterflyHawk = function(x) x$Family %in% c(bs,'Syrphidae',lh),
    BeesSyrphidsBeeFliesButterflyHawk = function(x) x$Family %in% combined
  )
}

## Return all records with usable and consistent labels, plus an audit of
## excluded specimens. Species-level and morphospecies labels both pass.
## Genus-level identifications are kept as genus-level nodes (never invented
## species). They should be treated carefully in downstream ecological metrics.
audit_network_ids <- function(spec.net) {
  required <- c('UniqueID','Site','Year','SampleRound', 'Family', 'Genus',
                'GenusSpecies', 'PlantGenusSpecies')
  if (length(missing <- setdiff(required, names(spec.net))))
    stop('Network data missing: ', paste(missing, collapse=', '))
  x <- as.data.frame(spec.net, stringsAsFactors=FALSE)
  for (nm in c('Site','Year','SampleRound','Family','Genus','GenusSpecies',
               'PlantGenusSpecies')) {
    x[[nm]] <- trimws(as.character(x[[nm]]))
    x[[nm]][is.na(x[[nm]])] <- ''
  }
  reason <- rep('', nrow(x))
  flag <- function(criterion, msg) {
    criterion[is.na(criterion)] <- FALSE
    reason[criterion] <<- ifelse(nzchar(reason[criterion]),
      paste(reason[criterion], msg, sep=';'), msg)
  }
  flag(!nzchar(x$Site) | !nzchar(x$Year) | !nzchar(x$SampleRound),
       'MISSING_EVENT_KEY')
  unknown_taxon <- function(z) tolower(z) %in% c(
    '', 'na', 'n/a', 'unknown', 'unidentified', 'unidentified insect',
    'unidentified plant', 'not identified', 'sp.', 'sp', 'spp.', 'spp',
    'insect', 'plant', 'none', '?')
  flag(unknown_taxon(x$Genus) | unknown_taxon(x$GenusSpecies),
       'MISSING_INSECT_ID')
  flag(unknown_taxon(x$PlantGenusSpecies), 'MISSING_PLANT_ID')
  flag(nzchar(x$Genus) & nzchar(x$GenusSpecies) &
    !(x$GenusSpecies == x$Genus |
      startsWith(x$GenusSpecies, paste0(x$Genus, ' '))),
       'INSECT_GENUS_LABEL_CONFLICT')

  ## Conflicting family assignments for the same genus should not be
  ## assigned to any guild until investigated.
  genus.fam <- unique(x[nzchar(x$Genus) & nzchar(x$Family),
                        c('Genus','Family'), drop=FALSE])
  fam.n <- table(genus.fam$Genus)
  conflicted <- names(fam.n[fam.n > 1])
  flag(x$Genus %in% conflicted, 'GENUS_FAMILY_CONFLICT')
  if ('InsectGBIFStatus' %in% names(x))
    flag(x$InsectGBIFStatus %in% 'REVIEW_FAMILY_CONFLICT',
         'GBIF_FAMILY_CONFLICT')
  ## A missing Family can be used in AllNetted; it cannot be inferred
  ## or silently assigned to a focal family based solely on a genus.

  ## Also flag plant genus contradictions when both names are available.
  if ('PlantGenus' %in% names(x)) {
    pg <- trimws(as.character(x$PlantGenus))
    pg[is.na(pg)] <- ''
    flag(nzchar(pg) & nzchar(x$PlantGenusSpecies) &
      !(x$PlantGenusSpecies == pg |
        startsWith(x$PlantGenusSpecies, paste0(pg, ' '))),
         'PLANT_GENUS_LABEL_CONFLICT')
  }
  x$NetworkIDReview <- reason
  list(usable=x[!nzchar(reason), ,drop=FALSE],
       rejected=x[nzchar(reason), ,drop=FALSE], all=x,
       conflicts=genus.fam[genus.fam$Genus %in% conflicted, ,drop=FALSE])
}

## Every network is plants (rows) x insects (columns).  Force plain
## numeric storage, even if the input is an 'xtabs' / 'table' object.
## Both dimnames and weighted interaction counts must be preserved.
plain_network_matrix <- function(web) {
  d <- dim(web)
  if (length(d) != 2L) stop('Expected a two-dimensional incidence table.')
  out <- matrix(as.numeric(web), nrow=d[[1L]], ncol=d[[2L]],
                dimnames=dimnames(web))
  if (any(!is.finite(out)) || any(out < 0))
    stop('Incidence matrix has missing, nonfinite or negative counts.')
  out
}

## Construct a proper dgCMatrix from positive cells ONLY, avoiding the
## ambiguous as(xtabs, 'dgCMatrix') / Matrix(xtabs) coercion path.
## This works even for a table whose 'xtabs' class has been retained by
## old or cached code, and avoids a dense->sparse class dispatch.
sparse_network_matrix <- function(web, binary=FALSE) {
  web <- plain_network_matrix(web)
  if (binary) web <- 1 * (web > 0)
  ij <- which(web > 0, arr.ind=TRUE)
  out <- Matrix::sparseMatrix(i=as.integer(ij[,1L]),
                              j=as.integer(ij[,2L]),
                              x=as.numeric(web[ij]),
                              dims=as.integer(dim(web)),
                              dimnames=dimnames(web), repr='C')
  if (!inherits(out, 'dgCMatrix'))
    stop('Internal error: sparse incidence is not dgCMatrix (class: ',
         paste(class(out), collapse=', '), ').')
  out
}

guild_incidence <- function(x) {
  if (!nrow(x)) return(matrix(0, nrow=0L, ncol=0L,
                              dimnames=list(character(),character())))
  web <- stats::xtabs(~PlantGenusSpecies + GenusSpecies, data=x,
                      drop.unused.levels=TRUE)
  plain_network_matrix(web)
}

## Keep graph conversion in one place, with an explicit sparse constructor.
## The binary graph has one edge per nonzero interaction, not one edge per
## individual specimen.
guild_network_graph <- function(web, weighted=TRUE) {
  if (!length(web)) return(igraph::make_empty_graph())
  incidence <- sparse_network_matrix(web, binary=!weighted)
  g <- igraph::graph_from_biadjacency_matrix(incidence,
              weighted=if (weighted) TRUE else NULL, directed=FALSE)
  igraph::V(g)$type <- ifelse(igraph::V(g)$type, 'Pollinator', 'Plant')
  g
}

make_guild_networks <- function(spec.net, event.keys, guilds, network_dir,
                                metrics_dir, audit_dir,
                                compute_bee_metrics=TRUE) {
  message('make_guild_networks: using xtabs-safe sparse conversion (v2)')
  dir.create(network_dir, recursive=TRUE, showWarnings=FALSE)
  dir.create(metrics_dir, recursive=TRUE, showWarnings=FALSE)
  dir.create(audit_dir, recursive=TRUE, showWarnings=FALSE)
  reviewed <- audit_network_ids(spec.net)
  utils::write.csv(reviewed$rejected,
      file.path(audit_dir, 'excluded_interactions.csv'), row.names=FALSE)
  utils::write.csv(reviewed$conflicts,
      file.path(audit_dir, 'genus_family_conflicts.csv'), row.names=FALSE)

  event.keys <- unique(as.data.frame(event.keys[,c('Site','Year','SampleRound'),
                                                 drop=FALSE]))
  for (nm in names(event.keys))
    event.keys[[nm]] <- trimws(as.character(event.keys[[nm]]))
  event.keys <- event.keys[stats::complete.cases(event.keys) &
     nzchar(event.keys$Site) & nzchar(event.keys$Year) &
     nzchar(event.keys$SampleRound), ,drop=FALSE]
  ## Include genuine observed events even if event.conditions lacks an
  ## entry. Empty networks are created only for verified event.keys.
  observed.keys <- unique(reviewed$usable[,c('Site','Year','SampleRound'),
                                           drop=FALSE])
  event.keys <- unique(rbind(event.keys, observed.keys))
  if (!nrow(event.keys)) stop('No identifiable network events.')
  event.keys <- event.keys[order(event.keys$Site,
                                 event.keys$Year,
                                 event.keys$SampleRound), ,drop=FALSE]
  event.keys$NetworkName <- paste(event.keys$Site, event.keys$Year,
                                  event.keys$SampleRound, sep='.')
  if (anyDuplicated(event.keys$NetworkName))
    stop('Network naming collision: Site/Year/SampleRound contain dots.')
  data <- reviewed$usable
  data$NetworkName <- paste(data$Site, data$Year,
                            data$SampleRound, sep='.')
  metrics.bees <- NULL
  all.manifests <- list()
  names.and.groups <- names(guilds)
  for (guild in names.and.groups) {
    matches <- guilds[[guild]](data)
    matches[is.na(matches)] <- FALSE
    guild.rows <- data[matches, ,drop=FALSE]
    per.event <- split(guild.rows, factor(guild.rows$NetworkName,
                            levels=event.keys$NetworkName), drop=FALSE)
    nets <- lapply(per.event, guild_incidence)
    names(nets) <- event.keys$NetworkName
    if (any(vapply(nets, function(z)
      any(!is.finite(z)) || any(z < 0), logical(1))))
      stop('Nonfinite network weights in ', guild)
    if (sum(vapply(nets, sum, numeric(1))) != nrow(guild.rows))
      stop('Network interactions do not reconcile with source specimens: ', guild)
    nets.graph <- lapply(nets, guild_network_graph, weighted=TRUE)
    nets.graph.uw <- lapply(nets, guild_network_graph, weighted=FALSE)
    sites <- event.keys$Site
    years <- event.keys$Year
    SRs <- event.keys$SampleRound
    names(sites) <- names(nets)
    names(years) <- names(nets)
    names(SRs) <- names(nets)
    manifest <- data.frame(Guild=guild,
      Site=sites, Year=years, SampleRound=SRs,
      NetworkName=names(nets),
      NInteractions=vapply(nets, sum, numeric(1)),
      NPlants=vapply(nets, nrow, integer(1)),
      NInsectTaxa=vapply(nets, ncol, integer(1)),
      stringsAsFactors=FALSE, row.names=NULL)
    all.manifests[[guild]] <- manifest
    utils::write.csv(manifest,
       file.path(audit_dir, paste0(guild, '_event_manifest.csv')),
       row.names=FALSE)
    save(nets, nets.graph, nets.graph.uw, sites, years, SRs, manifest,
      file=file.path(network_dir,
         sprintf('YearSR_PlantPollinator_%s.Rdata', guild)))
    if (guild == 'Bees' && compute_bee_metrics) {
      is.large <- vapply(nets, function(net) all(dim(net) > 2L), logical(1))
      if (any(is.large)) {
        ## Older site/round metric calculations only work on webs with
        ## at least three taxa on each side. Keep small matrices anyway.
        metrics.bees <- calcSpec(nets[is.large])
      } else {
        metrics.bees <- data.frame(Site=character(), Year=character(),
          SampleRound=character(), GenusSpecies=character(),
          speciesType=character(), SpSiteYear=character())
      }
      sp.lev <- metrics.bees
      save(sp.lev, file=file.path(metrics_dir,
           'YearSR_PlantPollinator_Bees.Rdata'))
    }
  }
  utils::write.csv(do.call(rbind, all.manifests),
    file.path(audit_dir, 'all_guild_event_manifest.csv'), row.names=FALSE)
  message('Produced ', length(guilds), ' named guild network lists for ',
          nrow(event.keys), ' sampling events (including documented zeros).')
  list(bee_metrics=metrics.bees,
       manifest=do.call(rbind, all.manifests),
       id_audit=reviewed)
}
