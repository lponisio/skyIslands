
check_unique_key <- function(dat, keys, label) {
  missing.cols <- setdiff(keys, names(dat))
  if (length(missing.cols)) stop(label, ': missing columns: ', paste(missing.cols, collapse=', '))
  nonempty <- stats::complete.cases(dat[, keys, drop=FALSE])
  for (k in keys) {
    if (is.character(dat[[k]]) || is.factor(dat[[k]]))
      nonempty <- nonempty & nzchar(trimws(as.character(dat[[k]])))
  }
  if (any(!nonempty)) {
    stop(label, ': ', sum(!nonempty), ' rows contain NA join keys (',
         paste(keys, collapse=', '), ').')
  }
  duplicates <- duplicated(dat[, keys, drop=FALSE]) |
    duplicated(dat[, keys, drop=FALSE], fromLast=TRUE)
  if (any(duplicates)) {
    print(utils::head(dat[duplicates, keys, drop=FALSE], 20))
    stop(label, ': duplicate keys for ', paste(keys, collapse=', '),
         ' (', sum(duplicates), ' rows). See preceding examples.')
  }
  invisible(TRUE)
}

check_rows_preserved <- function(before, after, label) {
  if (nrow(before) != nrow(after)) {
    stop(label, ': row count changed from ', nrow(before), ' to ', nrow(after),
         '; check join grain, duplicate keys, and missing keys.')
  }
  invisible(TRUE)
}
