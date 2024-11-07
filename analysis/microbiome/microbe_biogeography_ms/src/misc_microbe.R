
pdf.f <- function(f, file, ...) {
    cat(sprintf('Writing %s\n', file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
}

## add transparency to named colors
add.alpha <- function(col, alpha=0.2){
    apply(sapply(col, col2rgb)/255, 2,
          function(x)
              rgb(x[1], x[2], x[3],
                  alpha=alpha))
}


standardize <- function(x)
(x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)

unstandardize <- function(x, orig){
    (x*sd(orig, na.rm=TRUE)) + mean(orig, na.rm=TRUE)
}




standardize.axis <- function(x, orig)
(x-mean(orig, na.rm=TRUE))/sd(orig, na.rm=TRUE)


plot.res <- function(mod, mod.name){
    ## function tp plot diagnostic figures for mcmc
    pdf(sprintf("figures/diagnostics/%s_Diag.pdf", mod.name),
        height=11, width=8.5)
    plot(mod,  N = 4, ask = FALSE)
    dev.off()
}

## function to clean up white-space in a column of data (replaces all
## instances of white-space with " " and empty cells with ""
fix.white.space <- function(d) {
  d <- as.character(d)
  remove.first <- function(s) substr(s, 2, nchar(s))
  d <- gsub("      ", " ", d, fixed=TRUE)
  d <- gsub("     ", " ", d, fixed=TRUE)
  d <- gsub("    ", " ", d, fixed=TRUE)
  d <- gsub("   ", " ", d, fixed=TRUE)
  d <- gsub("  ", " ", d, fixed=TRUE)

  tmp <- strsplit(as.character(d), " ")
  d <- sapply(tmp, function(x) paste(x, collapse=" "))

  first <- substr(d, 1, 1)
  d[first==" "] <- remove.first(d[first==" "])
  d
}


## checks formula for any NAs
## NA check
check_for_NA <- function(which.form){
  for (x in which.form){
    anyNA <- length(spec.net[[x]][is.na(spec.net[[x]])])
    if (anyNA > 0){
    print(paste('microbe', x, "NAs:", length(spec.net[[x]][is.na(spec.net[[x]])])))
    #print(spec.net$GenusSpecies[spec.net[[x]][is.na(spec.net[[x]])]])
    }
  }
}

