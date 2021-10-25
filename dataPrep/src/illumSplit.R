library(parallel)

makeSpNet <- function(indiv.nets, sp.names, sp.abund){
    ## takes the indivual-level networks and takes the average across
    ## a species. then multiplies by the abundance at a site
    colnames(indiv.nets) <- sp.names
    int.mat <- sapply(unique(sp.names), function(i){
        if(dim(as.matrix(indiv.nets[,colnames(indiv.nets) == i]))[2] > 1){
            out <- apply(indiv.nets[,colnames(indiv.nets) == i], 1,
                         mean, na.rm=TRUE)
        } else {
            out <- indiv.nets[,colnames(indiv.nets) == i]
        }
    })

    abund <- sp.abund[order(match(colnames(int.mat),
                                  sp.abund$GenusSpecies)),]
    int.mat <- bipartite::empty(int.mat*abund$Abundance)

    return(int.mat)
}



illumSplit <- function(data,split,len){
    ## splits illumina data (which has already been attached
    ## to character data) by a split column, then converts the
    ## data into a list of network objects

    splitList <- lapply(unique(data[,split]), function(x){
        site <- subset(data,data[,split] == x)
        site.mx <- site[,len]
        site.mx <- apply(site.mx, 2, as.numeric)
        site.mx <- as.matrix(site.mx)
        if(ncol(site.mx) == 1){
            site.mx <- t(site.mx)

        }
        rownames(site.mx) <- site$UniqueID
        site.mx <- t(site.mx)
    })
    names(splitList) <- unique(data[,split])
    return(splitList)
}


netSums <- function(data, binary=FALSE, spec,
                    type="indiv", FUN=lapply, ...){
    ##takes a list of network objects (e.g., from illumSplit) and
    ##calculates summary stats using specieslevel. can be calculated
    ##either with reads acting as observations, or by converting the
    ##reads to presence/absence (binary = TRUE)

    ls <- FUN(names(data), function(x){
        data[[x]][is.na(data[[x]])] <- 0
        data[[x]] <- bipartite::empty(data[[x]])
        print(dim(data[[x]]))
        if(binary==FALSE){
            mets <- specieslevel(data[[x]])
            func.mets <- calcFuncUniqOrig(t(data[[x]]),
                                          rownames(data[[x]]),
                                          weights=rep(1,
                                                      length(rownames(data[[x]]))),
                                          ...)
            if(type == "indiv"){
                colnames(func.mets)[colnames(func.mets) == "GenusSpecies"]  <- "UniqueID"
            }
            ## df <- data.frame("UniqueID"=rownames(mets$`higher level`))
            df <- data.frame(mets$`higher level`)
            df$Site <- x
        } else {

            bin <- ifelse(data[[x]] > 0, 1, 0)

            mets <- specieslevel(bin)

            func.mets <- calcFuncUniqOrig(t(data[[x]]), rownames(data[[x]]),
                                          weights=rep(1,
                                          length(rownames(data[[x]]))),
                                          ...)

            df <- data.frame("UniqueID"=rownames(mets$`higher level`))
            df <- cbind(df,mets$`higher level`)
            df$Site <- x

        }
        if(type == "indiv"){
            df$GenusSpecies <- spec$GenusSpecies[match(rownames(df),
                                                       spec$UniqueID)]
            df$UniqueID <- rownames(df)
        } else {
            df$GenusSpecies <- rownames(df)
        }
        rownames(df) <- NULL
        df <- merge(df, func.mets)
        return(df)
    })

    do.call(rbind,ls)
}
