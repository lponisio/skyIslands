calcVar <- function(by.sp, var.fun){
    for.var <- by.sp[sapply(by.sp, nrow) > 1]
    met.var <- t(sapply(for.var, function(x) apply(x[,metrics], 2,
                                                 var.fun)))
    colnames(met.var) <- paste(colnames(met.var), "pv", sep=".")
    met.mean <- t(sapply(by.sp, function(x) apply(x[,metrics], 2,
                                                  mean, na.rm=TRUE)))
    colnames(met.mean) <- paste(colnames(met.mean), "mean", sep=".")
    species <- sapply(strsplit(rownames(met.mean), "[.]"), function(x){
        if(!is.na(x[3])) {
            out <- paste0(x[2], ".", x[3])
        }else{
            out <- x[2]
        }
        return(out)
    })
    year.site <- sapply(strsplit(rownames(met.mean), "[.]"),
                        function(x) x[1])
    var.met <- as.data.frame(met.mean)
    var.met <- cbind(var.met, met.var[match(rownames(var.met),
                                            rownames(met.var)),])
    rownames(var.met) <- NULL

    var.met$GenusSpecies <- species
    var.met$Year.Site <- year.site
    var.met$N <- sapply(by.sp, function(x) length(unique(x$Site)))
    maxN <- tapply(var.met$N, var.met$GenusSpecies, max)
    var.met$maxN <- maxN[match(var.met$GenusSpecies, names(maxN))]
    species.mets <- do.call(rbind, by.sp)
    var.met$speciesType <- species.mets$speciesType[match(var.met$GenusSpecies,
                                                    species.mets$GenusSpecies)]
    return(var.met)
}


PV <- function (Z) {
    n = length(Z)
    pairs = combn(Z, 2)
    min_z = apply(pairs, 2, min, na.rm=TRUE)
    max_z = apply(pairs, 2, max,  na.rm=TRUE)
    z = 1 - (min_z/max_z)
    PV = 2 * sum(z)/(n * (n - 1))
    return(PV)
}
