## Consistent taxonomy reconciliation for 1prep.R (insects, floral visits,
## quadrat vegetation, bloom surveys, and field-plant lookup).
## Explicitly pin the legacy GBIF BACKBONE: since rgbif 3.9 (Oct 2026)
## its unqualified name_backbone() uses COL Extended Release instead.
GBIF_BACKBONE_KEY <- 'd7dddbf4-2cf0-4f39-9b2a-bb099caae36c'

clean_taxon_string <- function(x) {
  x <- trimws(gsub('[[:space:]]+', ' ', as.character(x)))
  x[is.na(x) | x %in% c('NA', 'NaN', 'NULL', '""')] <- ''
  x
}

## Morphospecies, specimens identified to genus, and provisional IDs are
## explicitly preserved, NOT treated as spelling mistakes or queried as species.
is_provisional_taxon <- function(x) {
  x <- clean_taxon_string(x)
  tokens <- strsplit(x, ' ', fixed=TRUE)
  vapply(tokens, function(t) {
    if (length(t) < 2L) return(TRUE)
    epithet <- tolower(gsub('[.,]', '', t[-1]))
    any(grepl('^(sp|spp|morph|morpho|ms|cf|aff|nr|indet|undet|unknown|unidentified|group|complex|type|otu|asv)([0-9_-].*)?$', epithet)) ||
      any(grepl('^[0-9]+[a-z]?$', epithet)) ||
      any(grepl('^(sp|spp|morpho|morph)[._-]?[a-z0-9]+$', epithet))
  }, logical(1))
}

taxonomy_col <- function(x, col, default='') {
  if (is.data.frame(x) && col %in% names(x) && length(x[[col]]))
    return(as.character(x[[col]][1]))
  if (is.list(x) && !is.null(x[[col]]) && length(x[[col]]))
    return(as.character(x[[col]][1]))
  default
}

## A persisted, explicit correction audit: no data loss when GBIF is offline.
## Manual overrides are optional; only rows marked Approved == TRUE are applied.
reconcile_taxon_names <- function(names, kingdom, family=NULL,
                                  cache_file, overrides_file=NULL,
                                  refresh=FALSE) {
  stopifnot(length(names) == length(kingdom))
  if (is.null(family)) family <- rep('', length(names))
  names <- clean_taxon_string(names)
  kingdom <- clean_taxon_string(kingdom)
  family <- clean_taxon_string(family)
  x <- unique(data.frame(Kingdom=kingdom, OriginalName=names,
                         stringsAsFactors=FALSE))
  x <- x[nzchar(x$OriginalName), , drop=FALSE]
  dir.create(dirname(cache_file), recursive=TRUE, showWarnings=FALSE)
  cache.cols <- c('Kingdom','OriginalName','Backbone','MatchType','Confidence',
                  'Rank','CanonicalName','GBIFFamily','GBIFGenus','GBIFKingdom',
                  'UsageKey','LookupError')
  cache <- if (file.exists(cache_file))
    utils::read.csv(cache_file, stringsAsFactors=FALSE, na.strings=c(''))
    else as.data.frame(setNames(replicate(length(cache.cols), character(),
                                     simplify=FALSE), cache.cols))
  for (col in setdiff(cache.cols, names(cache))) cache[[col]] <- NA_character_
  cache <- cache[, cache.cols, drop=FALSE]
  cache$OriginalName <- clean_taxon_string(cache$OriginalName)
  cache$Kingdom <- clean_taxon_string(cache$Kingdom)
  cache <- cache[cache$Backbone == GBIF_BACKBONE_KEY &
                  !is.na(cache$Backbone), , drop=FALSE]
  cache <- cache[!duplicated(cache[c('Kingdom','OriginalName')]), ,drop=FALSE]

  x$Provisional <- is_provisional_taxon(x$OriginalName)
  query <- x[!x$Provisional, c('Kingdom','OriginalName'), drop=FALSE]
  lookup_key <- function(z) paste(z$Kingdom, z$OriginalName, sep='\r')
  prior <- match(lookup_key(query), lookup_key(cache))
  errored <- rep(FALSE, nrow(query))
  errored[!is.na(prior)] <- nzchar(clean_taxon_string(
    cache$LookupError[prior[!is.na(prior)]]))
  new <- query[is.na(prior) | errored | refresh, , drop=FALSE]
  can_lookup <- requireNamespace('rgbif', quietly=TRUE)
  if (!can_lookup && nrow(new)) warning('rgbif is not installed: ', nrow(new),
    ' names left unchanged and flagged UNCHECKED. Install rgbif for live matching.')
  if (can_lookup && nrow(new)) {
    message('GBIF Backbone: retrieving ', nrow(new), ' uncached taxonomic names.')
    consecutive.errors <- 0L
    for (i in seq_len(nrow(new))) {
      nm <- new$OriginalName[i]; kg <- new$Kingdom[i]
      err <- ''
      ans <- tryCatch(rgbif::name_backbone(name=nm, kingdom=kg,
                         checklistKey=GBIF_BACKBONE_KEY),
                      error=function(e) { err <<- conditionMessage(e); NULL })
      rec <- data.frame(Kingdom=kg, OriginalName=nm,
        Backbone=GBIF_BACKBONE_KEY,
        MatchType=taxonomy_col(ans,'matchType','NONE'),
        Confidence=taxonomy_col(ans,'confidence',''),
        Rank=taxonomy_col(ans,'rank',''),
        CanonicalName=taxonomy_col(ans,'canonicalName',''),
        GBIFFamily=taxonomy_col(ans,'family',''),
        GBIFGenus=taxonomy_col(ans,'genus',''),
        GBIFKingdom=taxonomy_col(ans,'kingdom',''),
        UsageKey=taxonomy_col(ans,'usageKey',''), LookupError=err,
        stringsAsFactors=FALSE)
      cache <- cache[lookup_key(cache) != lookup_key(rec), ,drop=FALSE]
      cache <- rbind(cache, rec)
      ## save incrementally: interrupted GBIF calls need not be repeated.
      utils::write.csv(cache, cache_file, row.names=FALSE, na='')
      if (nzchar(err)) {
        consecutive.errors <- consecutive.errors + 1L
        warning('GBIF lookup failed for ', kg, ' / ', nm, ': ', err)
        if (consecutive.errors >= 3L) {
          warning('GBIF unreachable: stopping lookups early; remaining names ',
                  'will stay unchanged and be marked UNCHECKED.')
          break
        }
      } else consecutive.errors <- 0L
    }
  }
  x <- merge(x, cache, by=c('Kingdom','OriginalName'), all.x=TRUE,
             sort=FALSE)
  x$Backbone <- GBIF_BACKBONE_KEY
  x$CorrectedName <- x$OriginalName
  x$Status <- ifelse(x$Provisional,'PROVISIONAL_PRESERVED','UNCHECKED')
  x$Confidence <- suppressWarnings(as.numeric(x$Confidence))
  x$MatchType <- clean_taxon_string(x$MatchType)
  x$Rank <- toupper(clean_taxon_string(x$Rank))
  x$CanonicalName <- clean_taxon_string(x$CanonicalName)
  x$GBIFFamily <- clean_taxon_string(x$GBIFFamily)
  x$GBIFGenus <- clean_taxon_string(x$GBIFGenus)
  x$GBIFKingdom <- clean_taxon_string(x$GBIFKingdom)
  names.parts <- strsplit(x$OriginalName, ' ', fixed=TRUE)
  x$OriginalGenus <- vapply(names.parts, `[`, character(1), 1)
  x$FamilyConflict <- FALSE
  ## A known family mismatch is a human-review item; do not automatically
  ## change either taxon or family (prevents false guild assignments).
  for (i in seq_len(nrow(x))) {
    this.families <- unique(family[kingdom == x$Kingdom[i] &
                            names == x$OriginalName[i]])
    this.families <- this.families[nzchar(this.families)]
    if (length(this.families) && nzchar(x$GBIFFamily[i]))
      x$FamilyConflict[i] <- any(tolower(this.families) !=
                                   tolower(x$GBIFFamily[i]))
  }
  found <- !x$Provisional & x$MatchType %in% c('EXACT','FUZZY') &
     x$Rank %in% c('SPECIES','SUBSPECIES','VARIETY','FORM','INFRASPECIFIC_NAME') &
     nzchar(x$CanonicalName) &
     (!nzchar(x$GBIFKingdom) | tolower(x$Kingdom) == tolower(x$GBIFKingdom))
  family_ok <- !x$FamilyConflict
  high_conf <- !is.na(x$Confidence) & x$Confidence >= 95
  ## Accepted-name synonyms may have a different genus. Only apply
  ## high-confidence, exact match synonyms automatically.
  exact <- found & x$MatchType == 'EXACT' & family_ok & high_conf
  ## For fuzzy matches additionally require the same genus and a small
  ## character edit distance: otherwise a real rare taxon may be rewritten.
  distances <- mapply(function(a,b) as.numeric(utils::adist(tolower(a),
    tolower(b))[1]), x$OriginalName, x$CanonicalName)
  fuzzy <- found & x$MatchType == 'FUZZY' & family_ok & high_conf &
    tolower(x$OriginalGenus) == tolower(x$GBIFGenus) & distances <= 2
  x$CorrectedName[exact | fuzzy] <- x$CanonicalName[exact | fuzzy]
  x$Status[exact] <- 'GBIF_EXACT_OR_SYNONYM'
  x$Status[fuzzy] <- 'GBIF_FUZZY_CORRECTED'
  unresolved <- !x$Provisional & !exact & !fuzzy
  x$Status[unresolved & !is.na(x$MatchType) & nzchar(x$MatchType)] <-
    'REVIEW_GBIF_MATCH'
  x$Status[x$FamilyConflict] <- 'REVIEW_FAMILY_CONFLICT'
  x$Status[!x$Provisional & !nzchar(x$MatchType)] <- 'UNCHECKED'

  ## Read the version-controlled approved corrections plus any additional
  ## investigator-approved corrections saved in the study checks directory.
  ## Both files use the same schema, and conflicting approvals are errors.
  if (!is.null(overrides_file)) {
    approved.paths <- overrides_file[file.exists(overrides_file)]
    if (length(approved.paths)) {
      read_approved <- function(p) {
        d <- utils::read.csv(p, stringsAsFactors=FALSE)
        needed <- c('Kingdom','OriginalName','Approved','CorrectedName')
        if (!all(needed %in% names(d)))
          stop('Taxon overrides in ', p, ' need: ',
               paste(needed, collapse=', '))
        d <- d[, needed, drop=FALSE]
        d$Kingdom <- clean_taxon_string(d$Kingdom)
        d$OriginalName <- clean_taxon_string(d$OriginalName)
        d$CorrectedName <- clean_taxon_string(d$CorrectedName)
        approved <- tolower(as.character(d$Approved)) %in%
          c('true','t','1','yes') & nzchar(d$CorrectedName)
        d[approved, ,drop=FALSE]
      }
      overrides <- do.call(rbind, lapply(approved.paths, read_approved))
      if (nrow(overrides)) {
        same.key <- paste(overrides$Kingdom, overrides$OriginalName, sep='\r')
        repeated <- duplicated(same.key) | duplicated(same.key, fromLast=TRUE)
        if (any(repeated)) {
          corrections <- split(overrides$CorrectedName[repeated], same.key[repeated])
          if (any(vapply(corrections, function(z) length(unique(z)) > 1L,
                         logical(1))))
            stop('Conflicting approved corrections in taxonomy override files.')
          overrides <- overrides[!duplicated(same.key), ,drop=FALSE]
        }
        idx <- match(lookup_key(x), lookup_key(overrides))
        selected <- !is.na(idx)
        x$CorrectedName[selected] <- overrides$CorrectedName[idx[selected]]
        x$Status[selected] <- 'MANUAL_APPROVED'
      }
    }
  }
  x
}

replace_from_taxonomy <- function(original, kingdom, lookup) {
  original <- clean_taxon_string(original)
  ind <- match(paste(kingdom, original, sep='\r'),
               paste(lookup$Kingdom, lookup$OriginalName, sep='\r'))
  updated <- lookup$CorrectedName[ind]
  updated[is.na(updated)] <- original[is.na(updated)]
  updated
}

## Optional, transparent study-level lumping retained for compatibility
## with earlier SI analyses. Never alter other plant genera.
si_plant_lump <- function(x, combine_erigeron=TRUE) {
  x <- clean_taxon_string(x)
  if (combine_erigeron)
    x[grepl('^Erigeron( |$)', x)] <- 'Erigeron spp.'
  x
}
