## Prepare specimen-level parasite assay indicators and richness.
##
## Input:  specimen data.frame from specimens-complete.csv, after taxon IDs.
## Output: list(spec = updated specimens, parasites = assay column names).
##
## Apidae == 1 means successful DNA extraction. Only those specimens can
## have documented parasite-negative (0) results assigned upstream in 1prep.R.
## Preserve explicit NA parasite calls from contaminated/indeterminate assays.
## For Apidae != 1 or NA, all parasite assay outcomes remain unknown (NA).
prepareParasites <- function(spec,
                             check_file = paste0(
                               '../../skyIslands_saved/data/checks/',
                               'screened_no_ID.csv')) {
  crithidias <- c('CrithidiaExpoeki', 'CrithidiaBombi',
                  'CrithidiaSpp', 'CrithidiaMellificae')
  nosemas <- c('NosemaBombi', 'NosemaCeranae')
  parasites <- c('AscosphaeraSpp', 'ApicystisSpp',
                 crithidias, nosemas)

  for (this.para in parasites) {
    spec[[this.para]] <- suppressWarnings(as.numeric(spec[[this.para]]))
  }
  spec$Apidae <- suppressWarnings(as.numeric(spec$Apidae))
  not.screened <- is.na(spec$Apidae) | spec$Apidae != 1
  spec[not.screened, parasites] <- NA
  spec$DNAExtractionSuccessful <- ifelse(spec$Apidae == 1, 1L, NA_integer_)

  ## A partially observed assay panel has indeterminate *total* richness.
  known.richness <- function(x) {
    if (anyNA(x)) return(NA_real_)
    sum(x)
  }

  ## Any positive test confirms presence; absence is confirmed only if
  ## every assay in the relevant group was observed and negative.
  known.presence <- function(x) {
    if (any(x > 0, na.rm=TRUE)) return(1L)
    if (anyNA(x)) return(NA_integer_)
    0L
  }

  spec$ParasiteRichness <- apply(spec[, parasites, drop=FALSE], 1,
                                 known.richness)
  spec$CrithidiaRichness <- apply(spec[, crithidias, drop=FALSE], 1,
                                  known.richness)
  spec$NosemaRichness <- apply(spec[, nosemas, drop=FALSE], 1,
                               known.richness)
  spec$PossibleParasite <- rowSums(!is.na(spec[, parasites, drop=FALSE]))
  spec$ParasitePresence <- apply(spec[, parasites, drop=FALSE], 1,
                                 known.presence)
  spec$CrithidiaPresence <- apply(spec[, crithidias, drop=FALSE], 1,
                                  known.presence)
  spec$NosemaPresence <- apply(spec[, nosemas, drop=FALSE], 1,
                               known.presence)

  ## Preserve original QC selection: a nonmissing Apidae result and
  ## a blank GenusSpecies. This includes non-success Apidae values, not
  ## exclusively Apidae == 1.
  check.spec <- spec[!is.na(spec$Apidae), , drop=FALSE]
  check.spec <- check.spec[check.spec$GenusSpecies == '', , drop=FALSE]
  dir.create(dirname(check_file), recursive=TRUE, showWarnings=FALSE)
  write.csv(check.spec, file=check_file, row.names=FALSE)

  ## Return the names too: the downstream parasite network uses `parasites`.
  list(spec=spec, parasites=parasites)
}
