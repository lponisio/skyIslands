
makeDataMultiLevel <- function(indiv.data, site.col, year.col="Year", weight.col.name="Weights",
                               site.id.col.name="SiteIDs"){
    ## split data by year
    indiv.data.split <- split(indiv.data, indiv.data[, year.col])

    ## maybe in the future this will need to be an sapply
    out.indiv.data <- lapply(indiv.data.split, addWeightCol,
                             site.col=site.col,  weight.col.name= weight.col.name,
                             site.id.col.name=site.id.col.name)

    out.indiv.data <- do.call(rbind, out.indiv.data)

    return(out.indiv.data)

}


addWeightCol <- function(each.year.dat, site.col, 
                         weight.col.name="Weights",
                         site.id.col.name="SiteIDs"){
    site.ids <- unlist(tapply(each.year.dat[, site.col],
                              each.year.dat[, site.col],
                              function(x) 1:length(x)))


    names(site.ids) <- NULL
    each.year.dat[, site.id.col.name] <- site.ids
    each.year.dat[, weight.col.name] <- each.year.dat[, site.id.col.name]
    each.year.dat[, weight.col.name][each.year.dat[weight.col.name] > 1] <- 0
    return(each.year.dat)
}
