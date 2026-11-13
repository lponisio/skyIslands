## Check the 2017 historical barcode crosswalk without freezing editable
## plant identifications. Called from relational/1prep.R, before assigning IDs.
##
## The crosswalk is ordered by historical TempID and maps original specimen
## records to physical barcode numbers. Changing FieldPlantID is legitimate
## manual taxonomic data cleaning, but inserting/reordering specimens is not.
validate_2017_barcode_map <- function(spec.data, bar.map, bbsl,
    expected_n_2017=20359L, expected_n_barcodes=9330L,
    audit_file='checks/taxonomy/2017_barcode_field_plant_edits.csv') {

  needed.spec <- c('Year','TempID','rawRowNum','Date','Site','SubSite',
                   'Collector','FieldPlantID','NetNumber','SampleRound',
                   'RawSpecimenCount')
  needed.cross <- c('TempID','BarcodeID','SourceRawRowNum','SourceDate',
                    'SourceSite','SourceSubSite','SourceCollector',
                    'SourceFieldPlantID','SourceNetNumber',
                    'SourceSampleRound','SourceSpecimenCount')
  if (!all(needed.spec %in% names(spec.data)))
    stop('Missing specimen barcode-audit columns: ',
         paste(setdiff(needed.spec,names(spec.data)),collapse=', '))
  if (!all(needed.cross %in% names(bar.map)))
    stop('Missing 2017 barcode crosswalk columns: ',
         paste(setdiff(needed.cross,names(bar.map)),collapse=', '))
  if (!'BarcodeID' %in% names(bbsl))
    stop('Missing BarcodeID in BBSL inventory.')

  rows.2017 <- which(as.character(spec.data$Year) == '2017')
  if (length(rows.2017) != expected_n_2017 ||
      nrow(bar.map) != expected_n_barcodes)
    stop('2017 specimen count or historical barcode count has changed. Review IDs first.')
  if (anyDuplicated(bar.map$TempID) || anyDuplicated(bar.map$BarcodeID) ||
      anyNA(bar.map$TempID) || anyNA(bar.map$BarcodeID) ||
      any(!nzchar(trimws(as.character(bar.map$BarcodeID)))))
    stop('The 2017 barcode reference contains duplicate or missing IDs.')

  first.2017 <- rows.2017[seq_len(nrow(bar.map))]
  if (!identical(as.integer(spec.data$TempID[first.2017]),
                 as.integer(bar.map$TempID)))
    stop('Historical 2017 TempID positions changed; do not shift barcodes.')

  normalize.audit.field <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x) | x %in% c('', 'NA')] <- ''
    x
  }

  ## Structural identity and collection provenance must still match exactly.
  ## The raw row number is especially important when collecting metadata
  ## happen to be identical among several specimens in the same net.
  immutable.fields <- c(SourceRawRowNum='rawRowNum', SourceDate='Date',
    SourceSite='Site', SourceSubSite='SubSite',
    SourceCollector='Collector', SourceNetNumber='NetNumber',
    SourceSampleRound='SampleRound', SourceSpecimenCount='RawSpecimenCount')
  for (field in names(immutable.fields)) {
    expected <- normalize.audit.field(bar.map[[field]])
    observed <- normalize.audit.field(
      spec.data[[immutable.fields[[field]]]][first.2017])
    bad <- which(expected != observed)
    if (length(bad))
      stop('2017 specimen order or collection metadata differs at ', field,
           ' (first affected TempID ', bar.map$TempID[bad[1]],
           '); barcode mapping intentionally not applied.')
  }

  if (nrow(bbsl) < nrow(bar.map) ||
      !identical(as.character(bbsl$BarcodeID[seq_len(nrow(bar.map))]),
                 as.character(bar.map$BarcodeID)))
    stop('BBSL inventory is not in historical order. Keep raw/BBSLSpecimens.csv unchanged.')

  ## Flower IDs can be corrected by hand after barcodes were assigned.
  ## Preserve the original value recorded at crosswalk creation in the audit.
  old.plant <- normalize.audit.field(bar.map$SourceFieldPlantID)
  new.plant <- normalize.audit.field(spec.data$FieldPlantID[first.2017])
  edited <- which(old.plant != new.plant)
  plant.edits <- data.frame(
    TempID=bar.map$TempID[edited],
    BarcodeID=bar.map$BarcodeID[edited],
    SourceRawRowNum=bar.map$SourceRawRowNum[edited],
    OriginalFieldPlantID=old.plant[edited],
    CurrentFieldPlantID=new.plant[edited],
    stringsAsFactors=FALSE)
  if (!is.null(audit_file)) {
    dir.create(dirname(audit_file), recursive=TRUE, showWarnings=FALSE)
    write.csv(plant.edits, audit_file, row.names=FALSE, na='')
  }
  if (length(edited))
    message('2017 barcode provenance verified; ',length(edited),
            ' historical specimen flower IDs have been edited. See ',audit_file)
  list(rows.2017=rows.2017, first.2017=first.2017,
       plant_edits=plant.edits)
}
